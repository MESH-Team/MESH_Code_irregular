NOT USED NOW
!subroutine mesh_diffusive(ppp,qqq, t0, t, tfin, saveInterval,j)
subroutine mesh_diffusive_forward(dtini_given, t0, t, tfin, saveInterval,j)


    use constants_module
    use arrays_module
    use matrix_module
    use var_module
    use arrays_section_module
    use xsec_attribute_module
    use subtools

    implicit none

    integer(kind=4), intent(in) :: j
    real(kind=4), intent(in) :: dtini_given, t0, t, tfin, saveInterval
    !doubleprecision, intent(in) :: t0, t, tfin, saveInterval


    real(kind=4) :: a1, a2, a3, a4, b1, b2, b3, b4, dd1, dd2, dd3, dd4, h1, h2, h3, h4, xt
    real(kind=4) :: qy, qxy, qxxy, qxxxy, ppi, qqi, rri, ssi, sxi, mannings, Sb, width

    real(kind=4) :: cour, cour2, q_sk_multi, sfi, r_interpol_time, r_interpo_nn, r_interpol, temp

    real(kind=4) :: y_norm_ds, y_crit_ds, S_ncomp, frds, area_0, width_0, hydR_0, errorY, currentQ
    integer :: tableLength

    integer(kind=4) :: i, pp


!!++++++++++++++++++++ Diffusive wave Forward sweep starts +++++++++++++++++++!!

        eei = -999.
		ffi = -999. !! What will be this value?
		exi = -999.
		fxi = -999.


		!!! steps for advection equation

        eei(1) = 1.0
		ffi(1) = 0. !! What will be this value?
		exi(1) = 0.
		fxi(1) = 0.

		ncomp = nx1(j)
            !print*, j, ncomp
            !pause
        do i = 2,ncomp


        !!!------ Calculation a1...a4, up to h4...
            cour = dtini / dx(i-1,j)
            cour2= abs( celerity(i,j) ) * cour

            a1 = 3.0 * cour2 ** 2.0 - 2.0 * cour2 ** 3.0
            a2 = 1 - a1
            a3 = ( cour2 ** 2.0 - cour2 ** 3.0 ) * dx(i-1,j)
            a4 = ( -1.0 * cour2 + 2.0 * cour2 ** 2.0 - cour2 ** 3.0 ) * dx(i-1,j)

            b1 = ( 6.0 * cour2 - 6.0 * cour2 ** 2.0 ) / ( -1.0 * dx(i-1,j) )
            b2 = - b1
            b3 = ( 2.0 * cour2 - 3.0 * cour2 ** 2.0 ) * ( -1.0 )
            b4 = ( -1.0 + 4.0 * cour2 - 3.0 * cour2 ** 2.0 ) * ( -1.0 )

            dd1 = ( 6.0 - 12.0 * cour2 ) / ( dx(i-1,j) ** 2.0 )
            dd2 = - dd1
            dd3 = ( 2.0 - 6.0 * cour2 ) / dx(i-1,j)
            dd4 = ( 4.0 - 6.0 * cour2 ) / dx(i-1,j)

            h1 = 12.0 / ( dx(i-1,j) ** 3.0 )
            h2 = - h1
            h3 = 6.0 / ( dx(i-1,j) ** 2.0 )
            h4 = h3


			qy   = a1 * oldQ(i-1,j) + a2 * oldQ(i,j) + a3 * qpx(i-1,j) + a4 * qpx(i,j)
			qxy  = b1 * oldQ(i-1,j) + b2 * oldQ(i,j) + b3 * qpx(i-1,j) + b4 * qpx(i,j)
			qxxy = dd1* oldQ(i-1,j) + dd2* oldQ(i,j) + dd3* qpx(i-1,j) + dd4* qpx(i,j)
			qxxxy= h1 * oldQ(i-1,j) + h2 * oldQ(i,j) + h3 * qpx(i-1,j) + h4 * qpx(i,j)

			!print*, qy, qxy, qxxy, qxxxy


			ppi = - theta * diffusivity(i,j) * dtini / ( dx(i-1,j) ** 2.0 )
			qqi = 1.0 - 2.0 * ppi
			rri = ppi
			ssi = qy  + dtini * diffusivity(i,j) * ( 1.0 - theta ) * qxxy + dtini * celerity(i,j) * lateralFlow(i)
			sxi = qxy + dtini * diffusivity(i,j) * ( 1.0 - theta ) * qxxxy+ dtini * celerity(i,j) * lateralFlow(i)/ dx(i-1,j)

			eei(i) = -1.0 * rri / ( ppi * eei(i-1) + qqi )                     !! copied from split operator method
			ffi(i) = ( ssi - ppi * ffi(i-1) ) / ( ppi * eei(i-1) + qqi )       !! copied from split operator method

			exi(i) = -1.0 * rri / ( ppi * exi(i-1) + qqi )
			fxi(i) = ( sxi - ppi * fxi(i-1) ) / ( ppi * exi(i-1) + qqi )

			!print*,i, dtini, celerity(i), lateralFlow(i,j)
        end do
        !pause 1002


			!pause 1001
        !! Applying d/s boundary
        ! qp(ncomp) = 0.
        qp(ncomp,j) = oldQ(ncomp-1,j)+lateralFlow(ncomp)*dx(ncomp-1,j)
        qpx(ncomp,j)= 0.

        do i = ncomp-1,1,-1

			qp(i,j) = eei(i) * qp(i+1,j) + ffi(i)
			qpx(i,j)= exi(i) *qpx(i+1,j) + fxi(i)

        end do
        qp(1,j) = newQ(1,j)
        newQ(1:ncomp,j) = qp(1:ncomp,j)
        dqp(1:ncomp,j) = newQ(1:ncomp,j)-oldQ(1:ncomp,j)
        dqc(1:ncomp,j) = dqp(1:ncomp,j)
        dap(1:ncomp,j) = 0.

       ! if (j .eq. 4) then
       !     print*, 'qp', qp(1:ncomp,j)
       !     print*, 'eei', eei(1:ncomp)
       !     print*, 'ffi', ffi(1:ncomp-1)
       !     print*, 'lateralFlow',lateralFlow(1:ncomp)
       !     print*, 'diffusivity',diffusivity(1:ncomp,j)
       !     print*, 'celerity',celerity(1:ncomp,j)
       ! end if
        !print*, '0', j, qp(1:ncomp,j)
        !print*, '1', j, celerity(1:ncomp,j)
        !print*, '2', j, diffusivity(1:ncomp,j)


!!++++++++++++++++++++ Diffusive wave Forward sweep ends +++++++++++++++++++!!
end subroutine mesh_diffusive_forward

subroutine mesh_diffusive_backward(dtini_given, t0, t, tfin, saveInterval,j)


    use constants_module
    use arrays_module
    use matrix_module
    use var_module
    use arrays_section_module
    use xsec_attribute_module
    use subtools

    implicit none

    integer(kind=4), intent(in) :: j
    real(kind=4), intent(in) :: dtini_given, t0, t, tfin, saveInterval
    !doubleprecision, intent(in) :: t0, t, tfin, saveInterval


    real(kind=4) :: a1, a2, a3, a4, b1, b2, b3, b4, dd1, dd2, dd3, dd4, h1, h2, h3, h4, xt
    real(kind=4) :: qy, qxy, qxxy, qxxxy, ppi, qqi, rri, ssi, sxi, mannings, Sb, width, slope

    real(kind=4) :: cour, cour2, q_sk_multi, sfi, r_interpol_time, r_interpo_nn, r_interpol, temp

    real(kind=4) :: y_norm_ds, y_crit_ds, S_ncomp, frds, area_0, width_0, hydR_0, errorY, currentQ
    integer :: tableLength, jj

    integer(kind=4) :: i, pp

!!++++++++++++++++++++ Diffusive wave Backward sweep starts +++++++++++++++++++!!

        !print*, j, 'normal depth', normalDepth(j)

        ncomp = nx1(j)

        S_ncomp = (-z(ncomp,j)+z(ncomp-1,j))/dx(ncomp-1,j)
        elevTable = xsec_tab(1,:,ncomp,j)
        areaTable = xsec_tab(2,:,ncomp,j)
        topwTable = xsec_tab(6,:,ncomp,j)

        newArea(ncomp,j) = r_interpol(elevTable,areaTable,nel,newY(ncomp,j))

        bo(ncomp,j) = r_interpol(elevTable,topwTable,nel,newY(ncomp,j))
        frds = abs(qp(ncomp,j))/sqrt(grav*( newArea(ncomp,j) )**3.0/bo(ncomp,j))



        !!! calculating y-normal and y-critical at downstream
        if (S_ncomp .gt. 0.) then
            currentQ = qp(ncomp,j)
            ! Calculating the Q_sk_multiplier for thie currentQ
            call calc_q_sk_multi(ncomp,j,currentQ,q_sk_multi)
            ! Calculating the normal depth and critical depth at the river reach upstream as an output
            call normal_crit_y(ncomp, j, q_sk_multi, S_ncomp, currentQ, y_norm_ds, y_crit_ds, area_norm, area_crit)
        end if

        !!! checking if d/s is critical
        !!! if d/s boundary is supercritical, the given boundary is ignored and normal depth is applied
        if (frds .ge. 1.0 ) newY(ncomp,j) = y_norm_ds


      !  if ((dimensionless_D(ncomp-1,j) .lt. 1) .and. (j .eq. nlinks)) then
      !      newY(ncomp,j) = y_norm_ds
      !  end if


        !if (j .eq. 8) print*,  j, ncomp,newY(ncomp,j)


        do i=ncomp,1,-1

            currentQ = qp(i,j)
            call calc_q_sk_multi(i,j,currentQ,q_sk_multi)

    !      Calculating : read all attributes from tab file
            elevTable = xsec_tab(1,:,i,j)
            convTable = xsec_tab(5,:,i,j)

            areaTable = xsec_tab(2,:,i,j)
            pereTable = xsec_tab(3,:,i,j)
            topwTable = xsec_tab(6,:,i,j)
    !     interpolate the cross section attributes based on water elevation
            xt=newY(i,j)

            currentCubicDepth=(elevTable-z(i,j))**3

            !co(i)  =q_sk_multi * r_interpol(elevTable,convTable,nel,xt)
            co(i)  =q_sk_multi * r_interpol(currentCubicDepth,convTable,nel,(newY(i,j)-z(i,j))**3.0)

            newArea(i,j) = r_interpol(elevTable,areaTable,nel,xt)
            pere(i,j) = r_interpol(elevTable,pereTable,nel,xt)
            bo(i,j) = r_interpol(elevTable,topwTable,nel,xt)

            sfi = ( qp(i,j) / co(i) ) ** 2.0


            ! Calculating empirical celerity:
            celerity2(i)=5.0 / 3.0 * sfi ** 0.3 * abs(qp(i,j)) ** 0.4 / bo(i,j) ** 0.4 / (1/(sk(i,j)*q_sk_multi)) ** 0.6
            diffusivity2(i) = abs(qp(i,j)) / 2.0 / bo(i,j) / sfi


!currentRoutingNormal(:,:), routingNotChanged(:,:)

            if (i .gt. 1) then
                !! If routing method is changed just a few time steps ago, we maintain the same routing to avoid oscillation
                if (routingNotChanged(i-1,j) .lt. minNotSwitchRouting) then
                    if (currentRoutingNormal(i-1,j) .eq. 1) then
                        slope = (z(i-1,j)-z(i,j))/dx(i-1,j)
                        call calc_q_sk_multi(i-1,j,qp(i-1,j),q_sk_multi)
                        ! applying normal depth to all the nodes
                        call normal_crit_y(i-1, j, q_sk_multi, slope, qp(i-1,j), newY(i-1,j), temp, newArea(i-1,j), temp)
                    else
                        newY(i-1,j) = newY(i,j) + sign ( sfi, qp(i,j) ) * dx(i-1,j)
                    end if
                else

                    !! If DSP: D is below 1.0, we switch to partial diffusive routing
                    if (dimensionless_D(i-1,j) .lt. 1) then
                        slope = (z(i-1,j)-z(i,j))/dx(i-1,j)
                        call calc_q_sk_multi(i-1,j,qp(i-1,j),q_sk_multi)
                        ! applying normal depth to all the nodes
                        call normal_crit_y(i-1, j, q_sk_multi, slope, qp(i-1,j), newY(i-1,j), temp, newArea(i-1,j), temp)

                         ! Book-keeping: changing from full diffusive to partial diffusive
                        if ( currentRoutingNormal(i-1,j) .eq. 0 ) routingNotChanged(i-1,j) = 0
                        currentRoutingNormal(i-1,j) = 1

                    !! If DSP: D is not below 1.0, we switch to full diffusive routing
                    else
                        newY(i-1,j) = newY(i,j) + sign ( sfi, qp(i,j) ) * dx(i-1,j)

                         ! Book-keeping: changing from partial diffusive to full diffusive
                        if ( currentRoutingNormal(i-1,j) .eq. 1 ) routingNotChanged(i-1,j) = 0
                        currentRoutingNormal(i-1,j) = 0

                    end if
                end if

                ! Book-keeping: Counting the number as for how many time steps the routing method is unchanged
                routingNotChanged(i-1,j) = routingNotChanged(i-1,j) + 1
            end if


       !if (j .eq. 4) print*, 'In diffusive backward 2', i


        end do

        !do i=1,ncomp
            !celerity(i) = sign ( sum(celerity2) / ncomp, qp(i) )

            celerity(1:ncomp,j) =  sum(celerity2(1:ncomp)) / ncomp

        !end do


		do i = 1, ncomp
            diffusivity(i,j)=sum(diffusivity2(1:ncomp)) / ncomp
			if (diffusivity(i,j) .gt. 10.) diffusivity(i,j) = 10. !!! Test
		end do



!!++++++++++++++++++++ Diffusive wave Backward sweep ends +++++++++++++++++++!!




end subroutine mesh_diffusive_backward
