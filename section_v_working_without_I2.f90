subroutine section(n)

    use constants_module
    use arrays_module
    use var_module
    use arrays_section_module

    implicit none

    ! Input
    integer, intent(in) :: n

    ! Locals
    integer :: i, pp, jj
    real(kind=4) :: beds, fs, hy, yyn, yyn_1, minZ, minZ_1, increm, increm_1
    real(kind=4) :: ar1, peri1, redi1, conv1, newI1, dpda_1, xt
    real(kind=4) :: eleNow, widthChange, avgEleNewSectionNow, r_interpol
    real(kind=4) :: waterDepthNow, cal_new_I2, avgWaterLevel

    do i=1,ncomp
        !print*, 'i=', i, 'n=', n
        depth(i)=y(n,i)-z(i)

!      Nazmul change: read all attributes from tab file

        write(file_num,'(i3.3)')i

        open(unit=19,file=trim(xSection_path)//file_num//'_tab')

        read(19,*)

        do pp=1,nel
          read(19,*,end=300)elevTable(pp),areaTable(pp),pereTable(pp),rediTable(pp),  &
                     convTable(pp),topwTable(pp),nwi1Table(pp),dPdATable(pp)
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
    end do

    do i=2,ncomp
        if(ityp(i-1) == 1) then

!---------------calculate I2 start--------------------------
      yyn=y(n,i)
      yyn_1=y(n,i-1)

      minZ=z(i)
      minZ_1=z(i-1)

      increm=(yyn-minZ)/(nel-1.0)
      increm_1=(yyn_1-minZ_1)/(nel-1.0)

      write(file_num_1,'(i3.3)')i-1
      xSection_tab_1=trim(xSection_path)//file_num_1//'_tab'

      open(unit=22,file=xSection_tab_1)

      read(22,*)

      do pp=1,maxTableLength
          read(22,*,end=500)elevTable_1(pp),ar1,peri1,redi1,conv1,  &
                  topwTable_1(pp),newI1, dpda_1
!          elevTable_1(pp)=el1
!          topwTable_1(pp)=tpW1
      enddo
500   close(22)

      do pp=1,nel
        eleNow=minZ+increm*(pp-1)
        topWNewSection(pp)=r_interpol(elevTable,topwTable,nel,elenow)
        eleNewSection(pp)=eleNow
        !print*, 'eleNewSection(pp)=',eleNewSection(pp)
        !pause 11

        eleNow=minZ_1+increm_1*(pp-1)
        topWNewSection_1(pp)=r_interpol(elevTable_1,topwTable_1,nel,elenow)
        eleNewSection_1(pp)=eleNow
        !print*, 'eleNewSection_1(pp)=',eleNewSection_1(pp)
        !pause 12
      end do

      cal_new_I2 = 0.0
      avgWaterLevel=(yyn+yyn_1)/2

      do pp=2,nel
        widthChange=topWNewSection(pp)-topWNewSection_1(pp)
        !print*, 'widthChange', widthChange
        avgEleNewSectionNow=(eleNewSection(pp)+eleNewSection_1(pp))/2
!       water depth at this point is the average depth of both sections
!       depth is taken up to the center of the slice
        waterDepthNow=avgWaterLevel-avgEleNewSectionNow+((increm+increm_1)/2)
        cal_new_I2=cal_new_I2+waterDepthNow*widthChange/dx(i)*((increm+increm_1)/2)
      enddo
      ci2(i)=cal_new_I2
      ci2(i)=0
            !ci2(i)=(depth(i)*depth(i)+depth(i-1)*depth(i-1))*(bo(i)-bo(i-1))/(dx(i-1)*4.)
!--------------calculate I2 end-----------------------------------

            beds=(z(i-1)-z(i))/dx(i-1)
            fs=f*0.5*q(n,i-1)*abs(q(n,i-1))/(co(i-1)*co(i-1))+f*0.5*q(n,i)*abs(q(n,i))/(co(i)*co(i))
            aso(i)=(area(i)+area(i-1))/2.0*(beds-fs)
            gso(i)=grav*(beds-fs)

            !Nazmul: TODO dbdx may needs to be changed
            dbdx(i)=(bo(i)-bo(i-1))/dx(i-1)
        end if
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
        !print*, xt, ' is not within the limit'
        !print*, 'maxval(x)= ', maxval(x), 'and minval(x)=', minval(x),'so',  xt, ' is not within the limit'
        if (xt.le. minval(x)) yt=minval(y)
        if (xt.ge. maxval(x)) yt=maxval(y)
    end if
    r_interpol = yt
    return
end function

