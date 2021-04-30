subroutine secpred(j)

    use constants_module
    use arrays_module
    use var_module
    use arrays_section_module
    use xsec_attribute_module
    use subtools

    implicit none

    integer, intent(in) :: j
    ! Locals
    integer :: i, pp, tableLength
    real :: beds, fs, hy, yyn, yyn_1, temp1, temp2, new_I2
    real :: xt, q_sk_multi, currentQ
    real :: r_interpo_nn

    !change 20191122
    ci1=0
    ci2=0


    do i=ncomp,1,-1
        if (applyNaturalSection .eq. 0) then
            yyn = areap(i,j) / bo(i,j) + z(i,j)
            depth(i) = areap(i,j) / bo(i,j)
            pere(i,j) = depth(i) * 2. + bo(i,j)
            hy = areap(i,j) / pere(i,j)
            co(i) = sk(i,j) * areap(i,j) * hy ** (2./3)
            ci1(i) = bo(i,j) * depth(i) *depth(i) / 2.
        else
        !print*, 'secpred', j, i
            if(i .lt. ncomp) then
                downstreamI2Tablep = I2Tablep
                downstreamSquareDepth = currentSquareDepth
                yyn_1 = yyn
            end if
! ----------------------------------------------------
            elevTable = xsec_tab(1,:,i,j)
            areaTable = xsec_tab(2,:,i,j)
            pereTable = xsec_tab(3,:,i,j)
            rediTable = xsec_tab(4,:,i,j)
            convTable = xsec_tab(5,:,i,j)
            topwTable = xsec_tab(6,:,i,j)
            nwi1Table = xsec_tab(7,:,i,j)
            dPdATable = xsec_tab(8,:,i,j)
            skkkTable = xsec_tab(11,:,i,j)
            currentSquareDepth=(elevTable-z(i,j))**2
            currentCubicDepth =(elevTable-z(i,j))**3
            I2Tablep = xsec_tab(9,:,i,j)
            I2Tablec = xsec_tab(10,:,i,j)

!      interpolate the cross section attributes based on predicted area
            xt=areap(i,j)

        ! NEW CHANGE: 20191119
        ! To get rid of negative predicted area
        !if (areap(i,j) .le. TOLERANCE) then
            !yyn = elevTable(1)+(elevTable(2)-elevTable(1))/100*5
            !areap(i) = 0.01
            !print*, 'At j = ', j, 'i=', i , ' pred area is negative'
        !else
        !print*, 'i=', i
        !print*, 'elevTable', elevTable
        !print*,  'areaTable', areaTable
        !print*, nel, xt,yyn
            call r_interpol(areaTable,elevTable,nel,xt,yyn)
            if (yyn .eq. -9999) then
                print*, 'At j = ',j,', i = ',i, 'interpolation of elevation was not possible'
                stop
            end if
        !end if

            depth(i) = yyn - z(i,j)
            xt = yyn

            call r_interpol(elevTable,pereTable,nel,xt,pere(i,j))
            if (pere(i,j) .eq. -9999) then
                print*, 'At j = ',j,', i = ',i, 'interpolation of peremeter was not possible'
                stop
            end if

            call r_interpol(elevTable,rediTable,nel,xt,hy)
            if (hy .eq. -9999) then
                print*, 'At j = ',j,', i = ',i, 'interpolation of hydraulic radius was not possible'
                stop
            end if

            call r_interpol(currentCubicDepth,convTable,nel,(xt-z(i,j))**3.0,co(i))
            if (co(i) .eq. -9999) then
                print*, 'At j = ',j,', i = ',i, 'interpolation of conveyance was not possible'
                stop
            end if

            call r_interpol(elevTable,topwTable,nel,xt,bo(i,j))
            if (bo(i,j) .eq. -9999) then
                print*, 'At j = ',j,', i = ',i, 'interpolation of top width was not possible'
                stop
            end if

            call r_interpol(currentSquareDepth,nwi1Table,nel,(depth(i))**2,ci1(i))
            if (ci1(i) .eq. -9999) then
                print*, 'At j = ',j,', i = ',i, 'interpolation of conveyance radius was not possible'
                stop
            end if

            call r_interpol(elevTable,skkkTable,nel,xt,sk(i,j))

        end if




!        write(*, 12) j,i,bo(i,j),depth(i),ci1(i)
!12      format(2i6, 3f24.17)



        !call r_interpol(elevTable,dPdATable,nel,xt,dpda(i))

        currentQ = qp(i,j)

        call calc_q_sk_multi(i,j,currentQ,q_sk_multi)

        co(i) = q_sk_multi*co(i)
! ----------------------------------------------------
        if(i .lt. ncomp) then
            if(ityp(i) == 1) then
            if (applyNaturalSection .eq. 0) then
                ci2(i)=(depth(i)*depth(i)+depth(i+1)*depth(i+1))*(bo(i+1,j)-bo(i,j))/(dx(i,j)*4.)
            else
! I2 opposite direction calculated as interpolation start
                xt=areap(i+1,j)
                call r_interpol(currentSquareDepth,I2Tablec,nel,(yyn-z(i,j))**2,temp1)
                call r_interpol(downstreamSquareDepth,downstreamI2Tablep, nel,(yyn_1-z(i+1,j))**2,temp2)
                new_I2 = (temp1+temp2)/2.0
                ci2(i)=new_I2
            end if
! I2 opposite direction calculated as interpolation end


                beds=(z(i,j)-z(i+1,j))/dx(i,j)
                fs=f*0.5*qp(i,j)*abs(qp(i,j))/(co(i)**2)+f*0.5*qp(i+1,j)*abs(qp(i+1,j))/(co(i+1)**2)
                aso(i,j)=(areap(i,j)+areap(i+1,j))/2.0*(beds-fs)
                gso(i,j)=grav*(beds-fs)
                !print*, j,i,beds,dx(i,j),aso(i,j),aso(i,j)
                dbdx(i,j)=(bo(i+1,j)-bo(i,j))/dx(i,j)
            end if
        end if

        !!! new for dkdh
   !     do pp = 2,nel
   !         if (yyn .le. elevTable(pp)) then
   !             dkdh(i)  =(convTable(pp)-convTable(pp-1))/(elevTable(pp)-elevTable(pp-1))
   !             EXIT
   !         endif
   !     end do
    end do
!pause
   ! print*, 'secpred'
   ! do i = 1, nx1(j)
   !     print*, i,depth(i),areap(i,j),ci1(i),co(i),ci2(i),aso(i,j),gso(i,j),bo(i,j)
   ! end do
   ! pause
    !gso(ncomp)=gso(ncomp-1)

end subroutine secpred
