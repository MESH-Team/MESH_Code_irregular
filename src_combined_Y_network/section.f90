subroutine section(j)

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


    !print*, 'section', j, oldY(1:ncomp,j),oldArea(1:ncomp,j)

    do i=1,ncomp


        depth(i)=oldY(i,j)-z(i,j)
        if (applyNaturalSection .eq. 0) then
            area(i) = depth(i) * bo(i,j)
            pere(i,j) = depth(i) * 2. + bo(i,j)
            hy = area(i) / pere(i,j)
            co(i) = sk(i,j) * area(i) * hy ** (2./3.)
            ci1(i) = bo (i,j) * depth(i) * depth(i) / 2.
        else
    !      Nazmul change: read all attributes from tab file
            if(i .gt. 1) then
                upstreamI2Tablec = I2Tablec
                upstreamSquareDepth = currentSquareDepth
            end if
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

        !     interpolate the cross section attributes based on water elevation
            xt=oldY(i,j)
            call r_interpol(elevTable,topwTable,nel,xt,bo(i,j))
            call r_interpol(elevTable,areaTable,nel,xt,area(i))
            call r_interpol(elevTable,pereTable,nel,xt,pere(i,j))
            call r_interpol(elevTable,rediTable,nel,xt,hy)
            call r_interpol(currentCubicDepth,convTable,nel,(xt-z(i,j))**3.0,co(i))
            call r_interpol(currentSquareDepth,nwi1Table,nel,(depth(i))**2.0,ci1(i))
            call r_interpol(elevTable,skkkTable,nel,xt,sk(i,j))
        !call r_interpol(elevTable,dPdATable,nel,xt,dpda(i))
        end if

        currentQ = oldQ(i,j)
        call calc_q_sk_multi(i,j,currentQ,q_sk_multi)
        co(i) = q_sk_multi*co(i)

! ----------------------------------------------------
        if(i .gt. 1) then
            if(ityp(i-1) == 1) then
            if (applyNaturalSection .eq. 0) then
                ci2(i)=(depth(i)*depth(i)+depth(i-1)*depth(i-1))*(bo(i,j)-bo(i-1,j))/(dx(i-1,j)*4.)
            else
! I2 calculated as interpolation start
                yyn=oldY(i,j)
                yyn_1=oldY(i-1,j)
                call r_interpol(currentSquareDepth,I2Tablep,nel,(depth(i))**2,temp1)
                call r_interpol(upstreamSquareDepth, upstreamI2Tablec, nel, (depth(i-1))**2,temp2)
                new_I2 = (temp1+temp2)/2.0
                ci2(i)=new_I2
! I2 calculated as interpolation end
            end if
                 beds=(z(i-1,j)-z(i,j))/dx(i-1,j)
                 fs=f*0.5*oldQ(i-1,j)*abs(oldQ(i-1,j))/(co(i-1)*co(i-1))+f*0.5*oldQ(i,j)*abs(oldQ(i,j))/(co(i)*co(i))
                 aso(i,j)=(area(i)+area(i-1))/2.0*(beds-fs)
                 gso(i,j)=grav*(beds-fs)
                 dbdx(i,j)=(bo(i,j)-bo(i-1,j))/dx(i-1,j)
            end if
        end if


        !!! new for dkdh
 !       do pp = 2,nel
 !           if (oldY(i) .le. elevTable(pp)) then
 !               dkdh(i)  =(convTable(pp)-convTable(pp-1))/(elevTable(pp)-elevTable(pp-1))
 !               EXIT
 !           endif
 !       end do
    end do
    !gso(1) = gso(2)
 !   print*, 'section'
 !   do i=1,ncomp
 !      print*, i,depth(i),area(i),ci1(i),co(i),ci2(i),aso(i,j),gso(i,j),bo(i,j)
!
!    end do
 !   pause


end subroutine section
