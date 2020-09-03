subroutine mesh_diffusive(ppp,qqq, t0, t, tfin, saveInterval)



    use constants_module
    use arrays_module
    use matrix_module
    use var_module
    use arrays_section_module
    use xsec_attribute_module
    use subtools

    implicit none

    integer(kind=4), intent(in) :: ppp, qqq
    real(kind=4), intent(in) :: t0, t, tfin, saveInterval
    !doubleprecision, intent(in) :: t0, t, tfin, saveInterval


    real(kind=4) :: a1, a2, a3, a4, b1, b2, b3, b4, dd1, dd2, dd3, dd4, h1, h2, h3, h4, xt
    real(kind=4) :: qy, qxy, qxxy, qxxxy, ppi, qqi, rri, ssi, sxi, mannings, Sb, width

    real(kind=4) :: cour, cour2, q_sk_multi, sfi, r_interpol_time, r_interpo_nn, r_interpol, temp, temp2, temp3

    real(kind=4) :: y_norm_ds, y_crit_ds, S_ncomp, frds, area_0, width_0, hydR_0, errorY, slope
    real :: area_n, errorX, normY
    integer :: tableLength

    integer(kind=4) :: i, pp

		call calculateDT(t0, t,saveInterval, cfl, tfin)

		lastKnownDiffuDT = dtini
        !print*, 'internal time', t, dtini
        !pause 5000

		!!! steps for advection equation

        eei(1) = 1.0
		ffi(1) = 0. !! What will be this value?
		exi(1) = 0.
		fxi(1) = 0.





        do i = 2,ncomp

                     !!!------ Calculation a1...a4, up to h4...
            cour = dtini / dx(i-1)
            cour2= abs( celerity(i) ) * cour

            a1 = 3.0 * cour2 ** 2.0 - 2.0 * cour2 ** 3.0
            a2 = 1 - a1
            a3 = ( cour2 ** 2.0 - cour2 ** 3.0 ) * dx(i-1)
            a4 = ( -1.0 * cour2 + 2.0 * cour2 ** 2.0 - cour2 ** 3.0 ) * dx(i-1)

            b1 = ( 6.0 * cour2 - 6.0 * cour2 ** 2.0 ) / ( -1.0 * dx(i-1) )
            b2 = - b1
            b3 = ( 2.0 * cour2 - 3.0 * cour2 ** 2.0 ) * ( -1.0 )
            b4 = ( -1.0 + 4.0 * cour2 - 3.0 * cour2 ** 2.0 ) * ( -1.0 )

            dd1 = ( 6.0 - 12.0 * cour2 ) / ( dx(i-1) ** 2.0 )
            dd2 = - dd1
            dd3 = ( 2.0 - 6.0 * cour2 ) / dx(i-1)
            dd4 = ( 4.0 - 6.0 * cour2 ) / dx(i-1)

            h1 = 12.0 / ( dx(i-1) ** 3.0 )
            h2 = - h1
            h3 = 6.0 / ( dx(i-1) ** 2.0 )
            h4 = h3


			qy   = a1 * oldQ(i-1) + a2 * oldQ(i) + a3 * qpx(i-1) + a4 * qpx(i)
			qxy  = b1 * oldQ(i-1) + b2 * oldQ(i) + b3 * qpx(i-1) + b4 * qpx(i)
			qxxy = dd1* oldQ(i-1) + dd2* oldQ(i) + dd3* qpx(i-1) + dd4* qpx(i)
			qxxxy= h1 * oldQ(i-1) + h2 * oldQ(i) + h3 * qpx(i-1) + h4 * qpx(i)

			!print*, qy, qxy, qxxy, qxxxy


			ppi = - theta * diffusivity(i) * dtini / ( dx(i-1) ** 2.0 )
			qqi = 1.0 - 2.0 * ppi
			rri = ppi
			ssi = qy  + dtini * diffusivity(i) * ( 1.0 - theta ) * qxxy + dtini * celerity(i) * lateralFlow(i)
			sxi = qxy + dtini * diffusivity(i) * ( 1.0 - theta ) * qxxxy+ dtini * celerity(i) * lateralFlow(i)/ dx(i-1)

			eei(i) = -1.0 * rri / ( ppi * eei(i-1) + qqi )                     !! copied from split operator method
			ffi(i) = ( ssi - ppi * ffi(i-1) ) / ( ppi * eei(i-1) + qqi )       !! copied from split operator method

			exi(i) = -1.0 * rri / ( ppi * exi(i-1) + qqi )
			fxi(i) = ( sxi - ppi * fxi(i-1) ) / ( ppi * exi(i-1) + qqi )

			!print*,i, dtini, celerity(i), lateralFlow(i)
        end do
        !pause 1002


			!pause 1001
        !! Applying d/s boundary
        ! qp(ncomp) = 0.
        qp(ncomp) = oldQ(ncomp-1)+lateralFlow(ncomp-1)*dx(ncomp-1)
        qpx(ncomp)= 0.

        do i = ncomp-1,1,-1

			qp(i) = eei(i) * qp(i+1) + ffi(i)
			qpx(i)= exi(i) *qpx(i+1) + fxi(i)
        end do

        S_ncomp = (-z(ncomp)+z(ncomp-1))/dx(ncomp-1)
        elevTable = xsec_tab(1,:,ncomp)
        areaTable = xsec_tab(2,:,ncomp)
        rediTable = xsec_tab(4,:,ncomp)
        topwTable = xsec_tab(6,:,ncomp)
        newArea(ncomp) = r_interpol(elevTable,areaTable,nel,newY(ncomp))
        bo(ncomp) = r_interpol(elevTable,topwTable,nel,newY(ncomp))
        frds = abs(qp(ncomp))/sqrt(grav*( newArea(ncomp) )**3.0/bo(ncomp))
        !print*, frds, newArea(ncomp), qp(ncomp), bo(ncomp)
        !pause 5002


        !!! calculating y-normal and y-critical at downstream
        if (S_ncomp .gt. 0.) then
              q_sk_multi = 1.0
              do pp = 1, size(Q_sk_tableEntry)
                if ( ( eachQSKtableNodeRange(1,pp) - ncomp) * ( eachQSKtableNodeRange(2,pp) - ncomp) .le. 0 ) then
                    tableLength = Q_sk_tableEntry(pp)
                    q_sk_multi = r_interpo_nn(Q_Sk_Table(1,1:tableLength,pp),Q_Sk_Table(2,1:tableLength,pp),tableLength,qp(ncomp))
                end if
              end do
              area_0 = newArea(ncomp)
              width_0= r_interpol(elevTable,topwTable,nel,newY(ncomp))
              area_crit=area_0
              errorY = 100.
            do while (errorY .gt. 0.0001)

              hydR_0 = r_interpol(areaTable,rediTable,nel,area_0)
              area_norm = qp(ncomp)/sk(ncomp)/q_sk_multi/ hydR_0 ** (2./3.) / sqrt(S_ncomp)

              errorY = abs(area_norm - area_0) / area_0
              area_0 = area_norm
              area_crit= (qp(ncomp) * qp(ncomp) * width_0 / grav) ** (1./3.)
              width_0  = r_interpol(areaTable,topwTable,nel,area_crit)
            enddo
          y_norm_ds = r_interpol(areaTable,elevTable,nel,area_0)
          y_crit_ds = r_interpol(areaTable,elevTable,nel,area_crit)
        end if
        !print*, y_norm_ds, y_crit_ds

        !!! checking if d/s is critical
        if (frds .ge. 1.0) newY(ncomp) = y_norm_ds

        qp(1) =r_interpol_time(USBoundary(1, 1:ppp),USBoundary(2, 1:ppp),ppp,t+dtini/60.)
        !newY(ncomp) =r_interpol_time(DSBoundary(1, 1:qqq),DSBoundary(2, 1:qqq),qqq,t+dtini/60.)
        depth(ncomp)=newY(ncomp)-z(ncomp)



        do i=ncomp,1,-1

            q_sk_multi = 1.0
            do pp = 1, size(Q_sk_tableEntry)
                if (  ( eachQSKtableNodeRange(1,pp) - i) * ( eachQSKtableNodeRange(2,pp) - i) .le. 0 ) then
                    q_sk_multi = r_interpo_nn(Q_Sk_Table(1,1:Q_sk_tableEntry(pp),pp),   &
                      Q_Sk_Table(2,1:Q_sk_tableEntry(pp),pp),Q_sk_tableEntry(pp),qp(i))
                end if
            end do

    !      Nazmul change: read all attributes from tab file
            elevTable = xsec_tab(1,:,i)
            convTable = xsec_tab(5,:,i)

            areaTable = xsec_tab(2,:,i)
            pereTable = xsec_tab(3,:,i)
            topwTable = xsec_tab(6,:,i)
    !     interpolate the cross section attributes based on water elevation
            xt=newY(i)
            co(i)  =q_sk_multi * r_interpol(elevTable,convTable,nel,xt)
            !width = r_interpol(elevTable,topwTable,nel,xt)



            newArea(i) = r_interpol(elevTable,areaTable,nel,xt)
            pere(i) = r_interpol(elevTable,pereTable,nel,xt)
            bo(i) = r_interpol(elevTable,topwTable,nel,xt)


            sfi = ( qp(i) / co(i) ) ** 2.0

            celerity2(i)=5.0 / 3.0 * sfi ** 0.3 * abs(qp(i)) ** 0.4 / bo(i) ** 0.4 / (1/(sk(i)*q_sk_multi)) ** 0.6
            diffusivity2(i) = abs(qp(i)) / 2.0 / bo(i) / sfi


        ! next 3 lines are commented out to check applicability of MC routing
        !    if (i .gt. 1) newY(i-1) = newY(i) + sign ( sfi, qp(i) ) * dx(i-1)
        !    if (i .eq. ncomp) print*,'qp=', qp
        !    if (i .gt. 1) print*, 'check', sfi, (z(i-1)-z(i))/dx(i-1), dx(i-1), newY(i)-z(i), newY(i-1)-z(i-1)




        ! for Muskingum Cunge routing:
        ! changing all WL as normal depth
            if (i .eq. ncomp) then
                slope = (z(i-1)-z(i))/dx(i-1)
            else
                slope = (z(i)-z(i+1))/dx(i)
            end if

            call normal_crit_y(i, q_sk_multi, slope, qp(i), newY(i), temp, newArea(i), temp)




        end do

        do i=1,ncomp
            !celerity(i) = sign ( sum(celerity2) / ncomp, qp(i) )

            celerity(i) =  sum(celerity2) / ncomp
            !celerity(i) =  celerity2(i) !!! Test
            !print*, celerity(i)
        end do


        !diffusivity = sum(diffusivity2) / ncomp

		do i = 1, ncomp
            diffusivity(i)=diffusivity2(i) !!! Test
			if (diffusivity(i) .gt. 10.) diffusivity(i) = 10. !!! Test
		end do
        ! Final update

        newQ = qp


        !diffusivity = 0.0


        !celerity = 0.5
        !diffusivity = 10.

        !pause 120


        !! Saving results
        !! print*, 'this time is complete'





end subroutine mesh_diffusive
