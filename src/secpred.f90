
subroutine secpred()

    use constants_module
    use arrays_module
    use var_module
    use arrays_section_module

     !  common/matrix/phi,theta,alfa2,alfa4,w
     !  common/sgate/ag,bg,option,yn,eps,dmeu,yw,bw,ydsn

    implicit none

    ! Locals
    integer :: i, pp, jj
    real(kind=4) :: beds, fs, hy, yyn, yyn_1, minZ, minZ_1, increm, increm_1
    real(kind=4) :: ar1, peri1, redi1, conv1, newI1, dpda_1, xt
    real(kind=4) :: eleNow, widthChange, avgEleNewSectionNow, r_interpol
    real(kind=4) :: waterDepthNow, cal_new_I2, avgWaterLevel


    do i=ncomp,1,-1
        !print*, 'secpred, i=', i
!        Nazmul: areap(i) is used through the section attribute table to get d(i)
        !if(i == 1) then
         !   print*, 'here'
        !end if

        if(i .lt. ncomp) then
            downstreamEleTable = elevTable
            downstreamTopwTable = topwTable
            downstreamAreaTable = areaTable
        end if
! ----------------------------------------------------
        write(file_num,'(i4.4)')i

        open(unit=19,file=trim(xSection_path)//file_num//'_tab')

        read(19,*)

        do pp=1,nel
            read(19,*,end=300)elevTable(pp),areaTable(pp),pereTable(pp),rediTable(pp),  &
                convTable(pp),topwTable(pp),nwi1Table(pp),dPdATable(pp)
        enddo
300     close(19)
        jj=pp-1

!      interpolate the cross section attributes based on predicted area
        xt=areap(i)

!      Nazmul: in this case, all the previous attributes are re-written
        pere(i)=r_interpol(areaTable,pereTable,nel,xt)
        hy     =r_interpol(areaTable,rediTable,nel,xt)
        co(i)  =r_interpol(areaTable,convTable,nel,xt)
        bo(i)  =r_interpol(areaTable,topwTable,nel,xt)
        ci1(i) =r_interpol(areaTable,nwi1Table,nel,xt)
        dpda(i)=r_interpol(areaTable,dPdATable,nel,xt)
        !print*, xt

        !pause 13
! ----------------------------------------------------

        !depth(i)=areap(i)/bo(i)
        !ci1(i)=bo(i)*depth(i)**2/2.0
        !hy=areap(i)/(2.0*depth(i)+bo(i))
        !co(i)=sk(i)*areap(i)*hy**(2.0/3.0)
        if(i .lt. ncomp) then
!---------------calculate I2 opposite direction start--------------------------
            yyn=r_interpol(areaTable,elevTable,nel,xt)
            xt=areap(i+1)
            yyn_1=r_interpol(downstreamAreaTable,downstreamEleTable,nel,xt)
            minZ=z(i)
            minZ_1=z(i+1)
            increm=(yyn-minZ)/(real(nel-1))
            increm_1=(yyn_1-minZ_1)/(real(nel-1))

            cal_new_I2 = 0
            avgWaterLevel=(yyn+yyn_1)/2

            do pp=2,nel
                eleNow=minZ+increm*(real(pp-1))
                topWNewSection(pp) = r_interpol(elevTable,topwTable,nel,eleNow)
                !pause 16
                eleNewSection(pp)=eleNow

                eleNow=minZ_1+increm_1*(real(pp-1))
                topWNewSection_downstream(pp) = r_interpol(downstreamEleTable,downstreamTopwTable,nel,eleNow)
                !pause 17
                eleNewSectionDownstream(pp)=eleNow

                widthChange=topWNewSection_downstream(pp)-topWNewSection(pp)
                avgEleNewSectionNow=(eleNewSection(pp)+eleNewSectionDownstream(pp))/2.0
        !       water depth at this point is the average depth of both sections
        !       depth is taken up to the center of the slice
                waterDepthNow=avgWaterLevel-avgEleNewSectionNow+((increm+increm_1)/4.0)

                cal_new_I2=cal_new_I2+waterDepthNow*widthChange/dx(i)*((increm+increm_1)/2.0)
            end do
!--------------calculate I2 opposite direction end-----------------------------------
            !print*,'at i=', i,'corr I2 =', cal_new_I2
            if(ityp(i) == 1) then
                ci2(i)=cal_new_I2
                beds=(z(i)-z(i+1))/dx(i)
                ! CHECK LOOP
                fs=f*0.5*qp(i)*abs(qp(i))/(co(i)**2)+f*0.5*qp(i+1)*abs(qp(i+1))/(co(i+1)**2)
                aso(i)=(areap(i)+areap(i+1))/2.0*(beds-fs)
                gso(i)=grav*(beds-fs)
                dbdx(i)=(bo(i+1)-bo(i))/dx(i)
            end if
        end if
    end do
    !pause 200

end subroutine secpred
