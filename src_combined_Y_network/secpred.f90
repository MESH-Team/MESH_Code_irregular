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
    real(kind=4) :: beds, fs, hy, yyn, yyn_1, temp1, temp2, new_I2
    real(kind=4) :: xt, q_sk_multi, currentQ
    real(kind=4) :: r_interpol, r_interpo_nn

    !change 20191122
    ci1=0
    ci2=0


    do i=ncomp,1,-1

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
        currentSquareDepth=(elevTable-z(i,j))**2
        I2Tablep = xsec_tab(9,:,i,j)
        I2Tablec = xsec_tab(10,:,i,j)

!      interpolate the cross section attributes based on predicted area
        xt=areap(i,j)

        ! NEW CHANGE: 20191119
        ! To get rid of negative predicted area
        if (areap(i,j) .le. TOLERANCE) then
            yyn = elevTable(1)+(elevTable(2)-elevTable(1))/100*5
            !areap(i) = 0.01
            print*, 'At i = ', i , ' pred area is negative'
        else
            yyn    =r_interpol(areaTable,elevTable,nel,xt)
        end if

        depth(i) = yyn - z(i,j)
        xt = yyn

        pere(i,j)=r_interpol(elevTable,pereTable,nel,xt)
        hy     =r_interpol(elevTable,rediTable,nel,xt)
        !co(i)  =r_interpol(elevTable,convTable,nel,xt)
        currentCubicDepth=(elevTable-z(i,j))**3
        co(i)  =q_sk_multi * r_interpol(currentCubicDepth,convTable,nel,(xt-z(i,j))**3.0)
        bo(i,j)  =r_interpol(elevTable,topwTable,nel,xt)
!        ci1(i) =r_interpol(areaTable,nwi1Table,nel,xt)
        ci1(i) =r_interpol(currentSquareDepth,nwi1Table,nel,(depth(i))**2)
        dpda(i)=r_interpol(elevTable,dPdATable,nel,xt)

        currentQ = qp(i,j)

!        q_sk_multi = 1.0
!        do pp = 1, noQSKtable(j)
!            if (  ( eachQSKtableNodeRange(1,pp,j) - i) * ( eachQSKtableNodeRange(2,pp,j) - i) .le. 0 ) then
!				tableLength = Q_sk_tableEntry(pp,j)
!                q_sk_multi = r_interpo_nn(Q_Sk_Table(1,1:tableLength,pp,j),Q_Sk_Table(2,1:tableLength,pp,j),tableLength,currentQ)
!            end if
!        end do

        call calc_q_sk_multi(i,j,currentQ,q_sk_multi)

        co(i) = q_sk_multi*co(i)
! ----------------------------------------------------
        if(i .lt. ncomp) then

! I2 opposite direction calculated as interpolation start
        xt=areap(i+1,j)
        temp1 = r_interpol(currentSquareDepth,I2Tablec,nel,(yyn-z(i,j))**2)
        temp2 = r_interpol(downstreamSquareDepth,downstreamI2Tablep, nel,(yyn_1-z(i+1,j))**2)
        new_I2 = (temp1+temp2)/2.0
! I2 opposite direction calculated as interpolation end

            if(ityp(i) == 1) then
                ci2(i)=new_I2
                beds=(z(i,j)-z(i+1,j))/dx(i,j)
                fs=f*0.5*qp(i,j)*abs(qp(i,j))/(co(i)**2)+f*0.5*qp(i+1,j)*abs(qp(i+1,j))/(co(i+1)**2)
                aso(i)=(areap(i,j)+areap(i+1,j))/2.0*(beds-fs)
                gso(i)=grav*(beds-fs)
                dbdx(i)=(bo(i+1,j)-bo(i,j))/dx(i,j)
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
    gso(ncomp)=gso(ncomp-1)

end subroutine secpred
