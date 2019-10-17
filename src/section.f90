subroutine section(n)

    use constants_module
    use arrays_module
    use var_module
    use arrays_section_module

    implicit none

    ! Input
    integer, intent(in) :: n

    ! Locals
    integer :: i, pp, jj, rez1
    real(kind=4) :: beds, fs, hy, yyn, yyn_1, minZ, minZ_1, increm, increm_1
    real(kind=4) :: ar1, peri1, redi1, conv1, newI1, dpda_1, xt
    real(kind=4) :: eleNow, widthChange, avgEleNewSectionNow, r_interpol
    real(kind=4) :: waterDepthNow, cal_new_I2, avgWaterLevel

    do i=1,ncomp

        if(i .gt. 1) then
            upstreamEleTable = elevTable
            upstreamTopwTable = topwTable
        end if

        !print*, 'i=', i, 'n=', n

        depth(i)=y(n,i)-z(i)

!      Nazmul change: read all attributes from tab file

        write(file_num,'(i3.3)')i

        open(unit=19,file=trim(xSection_path)//file_num//'_tab')

        read(19,*)

        do pp=1,nel
          read(19,*,end=300)elevTable(pp),areaTable(pp),pereTable(pp),rediTable(pp),  &
                     convTable(pp),topwTable(pp),nwi1Table(pp),dPdATable(pp)
          !print*, elevTable(pp),topwTable(pp)
        enddo
300     close(19)

        jj=pp-1
!     interpolate the cross section attributes based on water elevation
        xt=y(n,i)
!      bo = topwidth

        area(i)=r_interpol(elevTable,areaTable,nel,xt)
        pere(i)=r_interpol(elevTable,pereTable,nel,xt)
        hy     =r_interpol(elevTable,rediTable,nel,xt)
        co(i)  =r_interpol(elevTable,convTable,nel,xt)
        bo(i)  =r_interpol(elevTable,topwTable,nel,xt)
        ci1(i) =r_interpol(elevTable,nwi1Table,nel,xt)
        dpda(i)=r_interpol(elevTable,dPdATable,nel,xt)
        !pause 10

        !area(i)=bo(i)*depth(i)
        !ci1(i)=bo(i)*depth(i)**2/2.0
        !hy=area(i)/(2.0*depth(i)+bo(i))
        !co(i)=sk(i)*area(i)*hy**(2.0/3.0)
        ! print *, depth(i),area(i),ci1(i),y(n,i),z(i)
        if(i .gt. 1) then
    !---------------calculate I2 start--------------------------
            yyn=y(n,i)
            yyn_1=y(n,i-1)

            minZ=z(i)
            minZ_1=z(i-1)
            !print*, 'yyn', yyn, 'yyn_1', yyn_1, 'minZ', minZ, 'minZ_1', minZ_1

            increm=(yyn-minZ)/(nel-1.0)
            increm_1=(yyn_1-minZ_1)/(nel-1.0)
            !print*, 'increm', increm, 'increm_1', increm_1

            cal_new_I2 = 0.0
            avgWaterLevel=(yyn+yyn_1)/2

            do pp=2,nel
                !print*,'do I2, pp=',pp
                eleNow=minZ+increm*(pp-1.0)
                topWNewSection(pp)=r_interpol(elevTable,topwTable,nel,elenow)
                eleNewSection(pp)=eleNow
                !print*, 'eleNewSection(pp)=',eleNewSection(pp)
                !pause 11

                eleNow=minZ_1+increm_1*(pp-1)
                topWNewSection_upstream(pp)=r_interpol(upstreamEleTable,upstreamTopwTable,nel,elenow)
                eleNewSectionUpstream(pp)=eleNow
                !print*, 'eleNewSectionUpstream(pp)=',eleNewSectionUpstream(pp)
                !pause 12

                widthChange=topWNewSection(pp)-topWNewSection_upstream(pp)
                !print*, 'widthChange', widthChange
                avgEleNewSectionNow=(eleNewSection(pp)+eleNewSectionUpstream(pp))/2.0
                !       water depth at this point is the average depth of both sections
                !       depth is taken up to the center of the slice
                waterDepthNow=avgWaterLevel-avgEleNewSectionNow+((increm+increm_1)/4.0)
                cal_new_I2=cal_new_I2+waterDepthNow*widthChange/dx(i-1)*((increm+increm_1)/2.0)
                !print*, 'pp=', pp, 'I2=', cal_new_I2
            end do
    !--------------calculate I2 end-----------------------------------
            !print*,'at i=', i,'pred I2 =', cal_new_I2
            if(ityp(i-1) == 1) then
                 ci2(i)=cal_new_I2
                 beds=(z(i-1)-z(i))/dx(i-1)
                 fs=f*0.5*q(n,i-1)*abs(q(n,i-1))/(co(i-1)*co(i-1))+f*0.5*q(n,i)*abs(q(n,i))/(co(i)*co(i))
                 aso(i)=(area(i)+area(i-1))/2.0*(beds-fs)
                 gso(i)=grav*(beds-fs)
                !Nazmul: TODO dbdx may needs to be changed
                 dbdx(i)=(bo(i)-bo(i-1))/dx(i-1)
            end if
        end if
        !print*, 'section done'
    end do
    !pause 100

end subroutine section

function r_interpol(x,y,jj,xt)

    integer, intent(in) :: jj
    real, intent(in) :: xt, x(jj), y(jj)

    if (xt.lt.maxval(x) .and. xt.ge.minval(x)) then
        do j=1,jj-1
            if((x(j)-xt)*(x(j+1)-xt).le.0)then

               !print*,'j=',j,'j+1=',j+1
                !print*,'x(j)=',x(j),'x(j+1)=',x(j+1), 'xt=',xt
                yt=(xt-x(j))/(x(j+1)-x(j))*(y(j+1)-y(j))+y(j)
                !print*,'y(j)=',y(j),'y(j+1)=',y(j+1), 'yt=',yt
                EXIT
            endif
        end do
    else
        print*, xt, ' is not within the limit'
        print*, 'maxval(x)= ', maxval(x), 'and minval(x)=', minval(x),'so',  xt, ' is not within the limit'
        stop
        if (xt.le. minval(x)) yt=minval(y)
        if (xt.ge. maxval(x)) yt=maxval(y)
    end if
    r_interpol = yt
    return
end function

