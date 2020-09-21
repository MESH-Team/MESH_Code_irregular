subroutine mesh_diffusive(ppp,qqq, t0, t, tfin, saveInterval)



    use constants_module
    use arrays_module
    use matrix_module
    use var_module
    use arrays_section_module
    use xsec_attribute_module

    implicit none

    integer(kind=4), intent(in) :: ppp, qqq
    real(kind=4), intent(in) :: t0, t, tfin, saveInterval
    !doubleprecision, intent(in) :: t0, t, tfin, saveInterval


    real(kind=4) :: a1, a2, a3, a4, b1, b2, b3, b4, dd1, dd2, dd3, dd4, h1, h2, h3, h4, xt
    real(kind=4) :: qy, qxy, qxxy, qxxxy, ppi, qqi, rri, ssi, sxi, mannings, Sb, width

    real(kind=4) :: cour, cour2, q_sk_multi, sfi, r_interpol_time, r_interpo_nn, r_interpol, temp, temp2, temp3

    real(kind=4) :: y_norm_ds, y_crit_ds, S_ncomp, frds, area_0, width_0, hydR_0, errorY, slope, stg1, stg2
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
        qp(ncomp) = oldQ(ncomp-1)+lateralFlow(ncomp)*dx(ncomp-1)
        qpx(ncomp)= 0.

        do i = ncomp-1,1,-1
            !print*, qp(i), qpx(i), ffi(i)

			qp(i) = eei(i) * qp(i+1) + ffi(i)
			qpx(i)= exi(i) *qpx(i+1) + fxi(i)
        end do

        !print*, 'qp', qp

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
            call calc_q_sk_multi(ncomp,qp(ncomp),q_sk_multi)
            call normal_crit_y(ncomp, q_sk_multi, S_ncomp, qp(ncomp), y_norm_ds, y_crit_ds, area_norm, area_crit)
        end if
        !print*, y_norm_ds, y_crit_ds

        !!! checking if d/s is critical
        if (frds .ge. 1.0) newY(ncomp) = y_norm_ds

        qp(1) = newQ(1)
        !newY(ncomp) =r_interpol_time(DSBoundary(1, 1:qqq),DSBoundary(2, 1:qqq),qqq,t+dtini/60.)
        depth(ncomp)=newY(ncomp)-z(ncomp)

        !print*, 'check point 10'



        do i=ncomp,1,-1

            !q_sk_multi = 1.0
            !do pp = 1, size(Q_sk_tableEntry)
            !    if (  ( eachQSKtableNodeRange(1,pp) - i) * ( eachQSKtableNodeRange(2,pp) - i) .le. 0 ) then
            !        q_sk_multi = r_interpo_nn(Q_Sk_Table(1,1:Q_sk_tableEntry(pp),pp),   &
            !          Q_Sk_Table(2,1:Q_sk_tableEntry(pp),pp),Q_sk_tableEntry(pp),qp(i))
            !    end if
            !end do


            call calc_q_sk_multi(i,qp(i),q_sk_multi)







        ! for Muskingum Cunge routing:
        ! changing all WL as normal depth
            if (i .eq. ncomp) then
                slope = (z(i-1)-z(i))/dx(i-1)
            else
                slope = (z(i)-z(i+1))/dx(i)
            end if

            !print*, 'i = ', i, oldY(i), slope, qp(i)

            !call normal_crit_y(i, q_sk_multi, slope, qp(i), newY(i), temp, newArea(i), temp)

!      Nazmul change: read all attributes from tab file
            elevTable = xsec_tab(1,:,i)
            convTable = xsec_tab(5,:,i)

            areaTable = xsec_tab(2,:,i)
            pereTable = xsec_tab(3,:,i)
            topwTable = xsec_tab(6,:,i)
    !     interpolate the cross section attributes based on water elevation
            xt=newY(i)

            currentCubicDepth=(elevTable-z(i))**3

            !co(i)  = q_sk_multi * r_interpol(elevTable,convTable,nel,xt)
            co(i)  =q_sk_multi * r_interpol(currentCubicDepth,convTable,nel,(newY(i)-z(i))**3.0)

            bo(i) =  r_interpol(elevTable,topwTable,nel,xt)
            pere(i)= r_interpol(elevTable,pereTable,nel,xt)

            !width = r_interpol(elevTable,topwTable,nel,xt)



            sfi = ( qp(i) / co(i) ) ** 2.0


            celerity2(i)=5.0 / 3.0 * sfi ** 0.3 * abs(qp(i)) ** 0.4 / bo(i) ** 0.4 / (1/(sk(i)*q_sk_multi)) ** 0.6
            diffusivity2(i) = abs(qp(i)) / 2.0 / bo(i) / sfi

        if (i .gt. 1) then
			!! If DSP: D is below 1.0, we switch to partial diffusive routing
			if (dimensionless_D(i-1) .lt. 0.85) then
				slope = (z(i-1)-z(i))/dx(i-1)
				call calc_q_sk_multi(i-1,qp(i-1),q_sk_multi)
				! applying normal depth to all the nodes
				call normal_crit_y(i-1, q_sk_multi, slope, qp(i-1), newY(i-1), temp, newArea(i-1), temp)

				 ! Book-keeping: changing from full diffusive to partial diffusive
				if ( currentRoutingNormal(i-1) .ne. 1 ) routingNotChanged(i-1) = 0
				currentRoutingNormal(i-1) = 1

			!! If DSP: D is not below 1.0, we switch to full diffusive routing
			elseif ( (dimensionless_D(i-1) .ge. 0.85) .and. (dimensionless_D(i-1) .lt. 1.0) ) then
				slope = (z(i-1)-z(i))/dx(i-1)
				call calc_q_sk_multi(i-1,qp(i-1),q_sk_multi)
				! applying normal depth to all the nodes
				call normal_crit_y(i-1, q_sk_multi, slope, qp(i-1), stg1, temp, newArea(i-1), temp)

				stg2 = newY(i) + sign ( sfi, qp(i) ) * dx(i-1)
				newY(i-1) = ( stg2 * (dimensionless_D(i-1) - 0.85) + &
								stg1 * (1.0 - dimensionless_D(i-1)) ) / (1.0 - 0.85)
				 ! Book-keeping: changing from full diffusive to partial diffusive
				if ( currentRoutingNormal(i-1) .ne. 3 ) routingNotChanged(i-1) = 0
				currentRoutingNormal(i-1) = 3

			else
				newY(i-1) = newY(i) + sign ( sfi, qp(i) ) * dx(i-1)

				 ! Book-keeping: changing from partial diffusive to full diffusive
				if ( currentRoutingNormal(i-1) .ne. 0 ) routingNotChanged(i-1) = 0
				currentRoutingNormal(i-1) = 0

			end if
		end if

		! Book-keeping: Counting the number as for how many time steps the routing method is unchanged
		routingNotChanged(i-1) = routingNotChanged(i-1) + 1





        end do

        !do i=1,ncomp
            !celerity(i) = sign ( sum(celerity2) / ncomp, qp(i) )

            celerity(1:ncomp) =  sum(celerity2(1:ncomp)) / ncomp

        !end do


		do i = 1, ncomp
            diffusivity(i)=sum(diffusivity2(1:ncomp)) / ncomp
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

subroutine normal_crit_y(i, q_sk_multi, So, dsc, y_norm, y_crit, area_n, area_c)

    use constants_module
    use arrays_module
    use matrix_module
    use var_module
    use arrays_section_module
    use xsec_attribute_module

	implicit none

	integer, intent(in) :: i
	real, intent(in) :: q_sk_multi, So, dsc
	real, intent(out) :: y_norm, y_crit, area_n, area_c
	real :: area_0, width_0, errorY, hydR_0, r_interpol!, fro
	integer :: trapnm_app, recnm_app, iter


		elevTable = xsec_tab(1,:,i)
		areaTable = xsec_tab(2,:,i)
		rediTable = xsec_tab(4,:,i)
		topwTable = xsec_tab(6,:,i)
		area_0 = r_interpol(elevTable,areaTable,nel,oldY(i))
		width_0= r_interpol(elevTable,topwTable,nel,oldY(i))

		area_c=area_0
		errorY = 100.
		!pause
		do while (errorY .gt. 0.00001)

			hydR_0 = r_interpol(areaTable,rediTable,nel,area_0)
			area_n = dsc/sk(i)/q_sk_multi/ hydR_0 ** (2./3.) / sqrt(So)

			errorY = abs(area_n - area_0) / area_n
			area_0 = area_n
			area_c = (dsc * dsc * width_0 / grav) ** (1./3.)
			width_0  = r_interpol(areaTable,topwTable,nel,area_c)
			!fro=abs(dsc)/sqrt(grav*area_c**3.0/width_0)
		enddo

		y_norm = r_interpol(areaTable,elevTable,nel,area_0)
		y_crit = r_interpol(areaTable,elevTable,nel,area_c)


		!print*, 'check point 0', So, S_0,dsc,width_0, y_norm-z(i), y_crit-z(i), area_n, area_norm, area_c, fro
		!print*, 'check point -1', dsc, y_norm-z(1),  y_crit-z(1)
		!pause 500
end subroutine normal_crit_y


subroutine calc_q_sk_multi(i,currentQ,multipli)

    use arrays_module

	implicit none

	integer,intent(in) :: i
	real,intent(in) :: currentQ
	real,intent(out) :: multipli
	integer :: pp, tableLength
	real :: r_interpo_nn

	multipli = 1.0
	do pp = 1, size(Q_sk_tableEntry)
		if (  ( eachQSKtableNodeRange(1,pp) - i) * ( eachQSKtableNodeRange(2,pp) - i) .le. 0 ) then
			tableLength = Q_sk_tableEntry(pp)
			multipli = r_interpo_nn(Q_Sk_Table(1,1:tableLength,pp),Q_Sk_Table(2,1:tableLength,pp),tableLength,currentQ)
		end if
	end do

end subroutine calc_q_sk_multi

