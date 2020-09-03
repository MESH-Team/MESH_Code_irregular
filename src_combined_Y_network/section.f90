subroutine section(j)

    use constants_module
    use arrays_module
    use var_module
    use arrays_section_module
    use xsec_attribute_module

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

    do i=1,ncomp

        if(i .gt. 1) then
            upstreamI2Tablec = I2Tablec
            upstreamSquareDepth = currentSquareDepth
        end if

        depth(i)=oldY(i,j)-z(i,j)

!      Nazmul change: read all attributes from tab file
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

!     interpolate the cross section attributes based on water elevation
        xt=oldY(i,j)

        area(i)=r_interpol(elevTable,areaTable,nel,xt)
        !area(i)=r_interpol(x_tab,areaTable,nel,xt)
        pere(i,j)=r_interpol(elevTable,pereTable,nel,xt)
        hy     =r_interpol(elevTable,rediTable,nel,xt)
        co(i)  =r_interpol(elevTable,convTable,nel,xt)
        bo(i,j)  =r_interpol(elevTable,topwTable,nel,xt)
!        ci1(i) =r_interpol(elevTable,nwi1Table,nel,xt)
        ci1(i) =r_interpol(currentSquareDepth,nwi1Table,nel,(depth(i))**2)
        dpda(i)=r_interpol(elevTable,dPdATable,nel,xt)

        currentQ = oldQ(i,j)

        q_sk_multi = 1.0
        do pp = 1, noQSKtable(j)
            if (  ( eachQSKtableNodeRange(1,pp,j) - i) * ( eachQSKtableNodeRange(2,pp,j) - i) .le. 0 ) then
				tableLength = Q_sk_tableEntry(pp,j)
                q_sk_multi = r_interpo_nn(Q_Sk_Table(1,1:tableLength,pp,j),Q_Sk_Table(2,1:tableLength,pp,j),tableLength,currentQ)
            end if
        end do
        co(i) = q_sk_multi*co(i)

! ----------------------------------------------------
        if(i .gt. 1) then

! I2 calculated as interpolation start
        yyn=oldY(i,j)
        yyn_1=oldY(i-1,j)
        temp1 = r_interpol(currentSquareDepth,I2Tablep,nel,(depth(i))**2)
        temp2 = r_interpol(upstreamSquareDepth, upstreamI2Tablec, nel, (depth(i-1))**2)
        new_I2 = (temp1+temp2)/2.0
! I2 calculated as interpolation end

            if(ityp(i-1) == 1) then
                 ci2(i)=new_I2
                 beds=(z(i-1,j)-z(i,j))/dx(i-1,j)
                 fs=f*0.5*oldQ(i-1,j)*abs(oldQ(i-1,j))/(co(i-1)*co(i-1))+f*0.5*oldQ(i,j)*abs(oldQ(i,j))/(co(i)*co(i))
                 aso(i)=(area(i)+area(i-1))/2.0*(beds-fs)
                 gso(i)=grav*(beds-fs)
                 dbdx(i)=(bo(i,j)-bo(i-1,j))/dx(i-1,j)
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
    gso(1) = gso(2)




end subroutine section
