
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


    do i=1,ncomp
        print*, 'secpred, i=', i
!        Nazmul: areap(i) is used through the section attribute table to get d(i)
! ----------------------------------------------------
        write(file_num,'(i3.3)')i

        open(unit=19,file=trim(xSection_path)//file_num//'_tab')

        read(19,*)

        do pp=1,nel
            read(19,*,end=300)elevTable(pp),areaTable(pp),pereTable(pp),rediTable(pp),  &
                convTable(pp),topwTable(pp),nwi1Table(pp),dPdATable(pp)
!          elevTable(pp)=el1
!          areaTable(pp)=ar1
!          pereTable(pp)=peri1
!          rediTable(pp)=redi1
!          convTable(pp)=conv1
!          topwTable(pp)=tpW1
!          nwi1Table(pp)=newI1
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
    end do

    do i=1,ncomp-1
        if(ityp(i) == 1) then


!---------------calculate I2 opposite direction start--------------------------
      yyn=r_interpol(areaTable,elevTable,nel,xt)

      !pause 14

      write(file_num_1,'(i3.3)')i+1
      xSection_tab_1=trim(xSection_path)//file_num_1//'_tab'

      open(unit=21,file=xSection_tab_1)

      read(21,*)

      do pp=1,nel
          read(21,*,end=500)elevTable_1(pp),areaTable_1(pp),peri1,redi1,conv1,    &
                   topwTable_1(pp),newI1,dpda_1
      end do
500   close(21)

      xt=areap(i+1)
      yyn_1=r_interpol(areaTable_1,elevTable_1,nel,xt)
      !pause 15

      minZ=z(i)
      minZ_1=z(i+1)

!      ci2(i)=create_I2(xSection_path,i,yyn,yyn_1)
      increm=(yyn-minZ)/(real(nel-1))
      increm_1=(yyn_1-minZ_1)/(real(nel-1))

      do pp=1,nel
        eleNow=minZ+increm*(real(pp-1))
        topWNewSection(pp) = r_interpol(elevTable,topwTable,nel,eleNow)
        !pause 16
        eleNewSection(pp)=eleNow

        eleNow=minZ_1+increm_1*(real(pp-1))
        topWNewSection_1(pp) = r_interpol(elevTable_1,topwTable_1,nel,eleNow)
        !pause 17
        eleNewSection_1(pp)=eleNow
      end do

      cal_new_I2 = 0
      avgWaterLevel=(yyn+yyn_1)/2

      do pp=2,nel
        widthChange=topWNewSection_1(pp)-topWNewSection(pp)
        avgEleNewSectionNow=(eleNewSection(pp)+eleNewSection_1(pp))/2
!       water depth at this point is the average depth of both sections
!       depth is taken up to the center of the slice
        waterDepthNow=avgWaterLevel-avgEleNewSectionNow+((increm+increm_1)/2)

        cal_new_I2=cal_new_I2+waterDepthNow*widthChange/dx(i+1)*((increm+increm_1)/2)
      end do

!      IMPORTANT: here I2 is NOT calculated using the same previous formula, because
!      of the changes in sign and dx. The sign of widthchange is opposite and the
!      integral has to be devided by dx(i+1)
      ci2(i)=cal_new_I2
      ci2(i)=0
!--------------calculate I2 opposite direction end-----------------------------------

            !ci2(i)=(depth(i)*depth(i)+depth(i+1)*depth(i+1))*(bo(i+1)-bo(i))/(dx(i)*4.)
            beds=(z(i)-z(i+1))/dx(i)
            ! CHECK LOOP
            fs=f*0.5*qp(i)*abs(qp(i))/(co(i)**2)+f*0.5*qp(i+1)*abs(qp(i+1))/(co(i+1)**2)
            aso(i)=(areap(i)+areap(i+1))/2.0*(beds-fs)
            gso(i)=grav*(beds-fs)
          dbdx(i)=(bo(i+1)-bo(i))/dx(i)
        endif
    end do
    !pause 200

end subroutine secpred
