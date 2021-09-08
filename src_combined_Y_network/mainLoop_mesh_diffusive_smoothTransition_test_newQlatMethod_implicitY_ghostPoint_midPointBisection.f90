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

    integer, intent(in) :: j
    real, intent(in) :: dtini_given, t0, t, tfin, saveInterval
    !doubleprecision, intent(in) :: t0, t, tfin, saveInterval


    real :: a1, a2, a3, a4, b1, b2, b3, b4, dd1, dd2, dd3, dd4, h1, h2, h3, h4, xt, allqlat
    real :: qy, qxy, qxxy, qxxxy, ppi, qqi, rri, ssi, sxi, mannings, Sb, width

    real :: cour, cour2, q_sk_multi, sfi, r_interpol_time, r_interpo_nn, temp, alpha

    real :: y_norm_ds, y_crit_ds, S_ncomp, frds, area_0, width_0, hydR_0, errorY, currentQ
    integer :: tableLength, ll

    integer :: i, pp
    real :: eei_ghost, ffi_ghost, exi_ghost, fxi_ghost, qp_ghost, qpx_ghost


!!++++++++++++++++++++ Diffusive wave Forward sweep starts +++++++++++++++++++!!
!print*, 'dt=', dtini, 'cekerity=', celerity(1,1), 'diffusivity=', diffusivity(1,1)

! change 20210228: All qlat to a river reach is applied to the u/s boundary
! Note: lateralFlow(1,j) is already added to the boundary

        allqlat = sum(lateralFlow(2:ncomp,j) * dx(1:ncomp-1,j))

        !print*, 'j=',j
        !print*,lateralFlow(2:ncomp,j)
        !print*,dx(1:ncomp-1,j)
        !print*, 'allqlat',allqlat
        !pause 1010

        lateralFlow(:,j) = 0.


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

            if (i .eq. ncomp) then
                alpha = 1.0
            else
                alpha = dx(i,j) / dx(i-1,j)
            end if
            !alpha = 1.0

			qy   = a1 * oldQ(i-1,j) + a2 * oldQ(i,j) + a3 * qpx(i-1,j) + a4 * qpx(i,j)
			qxy  = b1 * oldQ(i-1,j) + b2 * oldQ(i,j) + b3 * qpx(i-1,j) + b4 * qpx(i,j)
			qxxy = dd1* oldQ(i-1,j) + dd2* oldQ(i,j) + dd3* qpx(i-1,j) + dd4* qpx(i,j)
			qxxxy= h1 * oldQ(i-1,j) + h2 * oldQ(i,j) + h3 * qpx(i-1,j) + h4 * qpx(i,j)

            ppi = - theta * diffusivity(i,j) * dtini / ( dx(i-1,j) ** 2.0 ) * 2.0 / (alpha*(alpha + 1.0)) * alpha
			qqi = 1.0 - ppi * (alpha + 1.0) / alpha
			rri = ppi / alpha

			ssi = qy  + dtini * diffusivity(i,j) * ( 1.0 - theta ) * qxxy !+ dtini * celerity(i,j) * lateralFlow(i,j)
			sxi = qxy + dtini * diffusivity(i,j) * ( 1.0 - theta ) * qxxxy !+ dtini * celerity(i,j) * lateralFlow(i)/ dx(i-1,j)

			eei(i) = -1.0 * rri / ( ppi * eei(i-1) + qqi )                     !! copied from split operator method
			ffi(i) = ( ssi - ppi * ffi(i-1) ) / ( ppi * eei(i-1) + qqi )       !! copied from split operator method

			exi(i) = -1.0 * rri / ( ppi * exi(i-1) + qqi )
			fxi(i) = ( sxi - ppi * fxi(i-1) ) / ( ppi * exi(i-1) + qqi )

        end do

        !!! Ghost point calculation start

            cour = dtini / dx(ncomp-1,j)
            cour2= abs( celerity(ncomp-1,j) ) * cour

            a1 = 3.0 * cour2 ** 2.0 - 2.0 * cour2 ** 3.0
            a2 = 1 - a1
            a3 = ( cour2 ** 2.0 - cour2 ** 3.0 ) * dx(ncomp-1,j)
            a4 = ( -1.0 * cour2 + 2.0 * cour2 ** 2.0 - cour2 ** 3.0 ) * dx(ncomp-1,j)

            b1 = ( 6.0 * cour2 - 6.0 * cour2 ** 2.0 ) / ( -1.0 * dx(ncomp-1,j) )
            b2 = - b1
            b3 = ( 2.0 * cour2 - 3.0 * cour2 ** 2.0 ) * ( -1.0 )
            b4 = ( -1.0 + 4.0 * cour2 - 3.0 * cour2 ** 2.0 ) * ( -1.0 )

            dd1 = ( 6.0 - 12.0 * cour2 ) / ( dx(ncomp-1,j) ** 2.0 )
            dd2 = - dd1
            dd3 = ( 2.0 - 6.0 * cour2 ) / dx(ncomp-1,j)
            dd4 = ( 4.0 - 6.0 * cour2 ) / dx(ncomp-1,j)

            h1 = 12.0 / ( dx(ncomp-1,j) ** 3.0 )
            h2 = - h1
            h3 = 6.0 / ( dx(ncomp-1,j) ** 2.0 )
            h4 = h3

            alpha = 1.0

			qy   = a1 * oldQ(ncomp,j) + a2 * oldQ(ncomp-1,j) + a3 * qpx(ncomp,j) + a4 * qpx(ncomp-1,j)
			qxy  = b1 * oldQ(ncomp,j) + b2 * oldQ(ncomp-1,j) + b3 * qpx(ncomp,j) + b4 * qpx(ncomp-1,j)
			qxxy = dd1* oldQ(ncomp,j) + dd2* oldQ(ncomp-1,j) + dd3* qpx(ncomp,j) + dd4* qpx(ncomp-1,j)
			qxxxy= h1 * oldQ(ncomp,j) + h2 * oldQ(ncomp-1,j) + h3 * qpx(ncomp,j) + h4 * qpx(ncomp-1,j)
			!print*, qy, qxy, qxxy, qxxxy

            ppi = - theta * diffusivity(ncomp,j) * dtini / ( dx(ncomp-1,j) ** 2.0 ) * 2.0 / (alpha*(alpha + 1.0)) * alpha
			qqi = 1.0 - ppi * (alpha + 1.0) / alpha
			rri = ppi / alpha

			ssi = qy  + dtini * diffusivity(ncomp-1,j) * ( 1.0 - theta ) * qxxy !+ dtini * celerity(i,j) * lateralFlow(i,j)
			sxi = qxy + dtini * diffusivity(ncomp-1,j) * ( 1.0 - theta ) * qxxxy !+ dtini * celerity(i,j) * lateralFlow(i)/ dx(i-1,j)


			eei_ghost = -1.0 * rri / ( ppi * eei(ncomp) + qqi )                     !! copied from split operator method
			ffi_ghost = ( ssi - ppi * ffi(ncomp) ) / ( ppi * eei(ncomp) + qqi )       !! copied from split operator method

			exi_ghost = -1.0 * rri / ( ppi * exi(ncomp) + qqi )
			fxi_ghost = ( sxi - ppi * fxi(ncomp) ) / ( ppi * exi(ncomp) + qqi )

			!print*, eei(ncomp),eei_ghost, ffi(ncomp), ffi_ghost

        !!! Ghost point calculation end


    !pause 1001
        !! Applying d/s boundary
        ! qp(ncomp) = 0.
        !qp(ncomp,j) = oldQ(ncomp-1,j)+lateralFlow(ncomp,j)*dx(ncomp-1,j)
        !qpx(ncomp,j)= 0.

        qp_ghost = oldQ(ncomp-1,j)
        qpx_ghost= 0.

        qp(ncomp,j) = eei(ncomp) * qp_ghost + ffi(ncomp)
        qpx(ncomp,j)= exi(ncomp) *qpx_ghost + fxi(ncomp)
        !print*, eei_ghost, ffi_ghost; pause

        do i = ncomp-1,1,-1

			qp(i,j) = eei(i) * qp(i+1,j) + ffi(i)
			qpx(i,j)= exi(i) *qpx(i+1,j) + fxi(i)

			!print*, i, qp(i,j), qpx(i,j)

        end do
        !qp(1:9,1)=1.0; qp(10:ncomp,1)=11.0;
        !pause
        !if (j .eq. 264) print*, 'qp',(qp(i,j),  i=1, ncomp),'ql', (lateralFlow(i,j),  i=1, ncomp),'oldq', (oldQ(i,j), i=1,ncomp)
        !if (j .eq. 264) pause 6500
        qp(1,j) = newQ(1,j)

        !if (j .eq. 277) then
        !    print*, 'qp2', t, qp(1,j)
        !end if

        ! change 20210228: All qlat to a river reach is applied to the u/s boundary
        qp(1,j) = qp(1,j) + allqlat

        !if (j .eq. 277) then
        !    print*, 'qp3', t, qp(1,j)
        !end if

        do i=1,ncomp
            if (abs(qp(i,j)) .lt. 0.02831) then
                qp(i,j) = 0.02831
            end if
        end do

        newQ(1:ncomp,j) = qp(1:ncomp,j)
        dqp(1:ncomp,j) = newQ(1:ncomp,j)-oldQ(1:ncomp,j)
        dqc(1:ncomp,j) = dqp(1:ncomp,j)
        dap(1:ncomp,j) = 0.
!do i=1,ncomp
!    print*, i, qp(i,j), qpx(i,j)
!end do
!pause

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

    integer, intent(in) :: j
    real, intent(in) :: dtini_given, t0, t, tfin, saveInterval
    !doubleprecision, intent(in) :: t0, t, tfin, saveInterval


    real :: a1, a2, a3, a4, b1, b2, b3, b4, dd1, dd2, dd3, dd4, h1, h2, h3, h4, xt
    real :: qy, qxy, qxxy, qxxxy, ppi, qqi, rri, ssi, sxi, mannings, Sb, width, slope

    real :: cour, cour2, q_sk_multi, sfi, r_interpol_time, r_interpo_nn, temp, dkdh
    real :: D_lim1, D_lim2, y_norm, y_crit, area_n, area_c, chnWidth, vel

    real :: y_norm_ds, y_crit_ds, S_ncomp, frds, area_0, width_0, hydR_0, errorY, currentQ, stg1, stg2
    integer :: tableLength, jj, newMassBalance, iii


    real :: elevTable_1(nel),areaTable_1(nel),rediTable_1(nel),convTable_1(nel),topwTable_1(nel),currentSquaredDepth_1(nel)
    real :: skkkTable_1(nel)
    real :: pereTable_1(nel),depthYi,tempDepthi_1,tempCo_1,temp_q_sk_multi_1,tempY_1,tempArea_1,tempRadi_1,tempbo_1,ffy
    real :: ffy_1, ffy_2, ffy_3, tempCo_2, tempCo_3, tempsfi_2, tempsfi_3
    real :: ffprime,tempDepthi_1_new,tempsfi_1,toll, dkda, tempPere_1, tempsk_1, tempY_2, tempY_3, usFroud
    integer :: depthCalOk(ncomp), newtonRaphson

    integer :: i, pp

!!++++++++++++++++++++ Diffusive wave Backward sweep starts +++++++++++++++++++!!

        !print*, j, 'normal depth', normalDepth(j)
        D_lim1 = -10.
        D_lim2 = -15.
        newtonRaphson = 1

        ncomp = nx1(j)


        S_ncomp = (-z(ncomp,j)+z(ncomp-1,j))/dx(ncomp-1,j)
        elevTable = xsec_tab(1,:,ncomp,j)
        areaTable = xsec_tab(2,:,ncomp,j)
        topwTable = xsec_tab(6,:,ncomp,j)

        depthCalOk(ncomp) = 1



        call r_interpol(elevTable,areaTable,nel,newY(ncomp,j),newArea(ncomp,j))
		if (newArea(ncomp,j) .eq. -9999) then
			print*, 'At j = ',j,', i = ',ncomp, 'time =',t, 'interpolation of newArea was not possible'
			stop
		end if

        call r_interpol(elevTable,topwTable,nel,newY(ncomp,j),bo(ncomp,j))
		if (bo(ncomp,j) .eq. -9999) then
			print*, 'At j = ',j,', i = ',ncomp, 'time =',t, 'interpolation of bo was not possible'
			stop
		end if


        !!! calculating y-normal and y-critical at downstream
        !if (j .ne. 2) then
 !       if (S_ncomp * qp(ncomp,j) .gt. 0.) then
 !           currentQ = qp(ncomp,j)
            ! Calculating the Q_sk_multiplier for thie currentQ
 !           call calc_q_sk_multi(ncomp,j,currentQ,q_sk_multi)
            !print*, '102', j, q_sk_multi
            ! Calculating the normal depth and critical depth at the river reach upstream as an output
 !           call normal_crit_y(ncomp, j, q_sk_multi, S_ncomp, currentQ, y_norm_ds, y_crit_ds, area_norm, area_crit)

  !          frds = abs(qp(ncomp,j))/sqrt(grav*( newArea(ncomp,j) )**3.0/bo(ncomp,j))
            !!! checking if d/s is critical
            !!! if d/s boundary is supercritical, the given boundary is ignored and normal depth is applied
  !          if (frds .ge. 1.0 ) newY(ncomp,j) = y_norm_ds
  !          if (frds .ge. 1.0 ) newArea(ncomp,j) = area_norm


            !!! temp !!!
  !          normalDepthAtNodes(ncomp,j) = y_norm_ds
            !!! temp !!!


  !      end if
        !end if

		!print*, '1001',ncomp, j,frds, q_sk_multi, S_ncomp, currentQ, y_norm_ds, y_crit_ds, area_norm, area_crit


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
            skkkTable = xsec_tab(11,:,i,j)
    !     interpolate the cross section attributes based on water elevation
            xt=newY(i,j)

            currentSquareDepth=(elevTable-z(i,j))**2.


            !pause 5000
            call r_interpol(currentSquareDepth,convTable,nel,(newY(i,j)-z(i,j))**2.0,co(i))
            !if (j .eq. 3) print*, i,'(newY(i,j)-z(i,j))**3.0',(newY(i,j)-z(i,j))**3.0, 'co', co(i)
            !pause
            if (co(i) .eq. -9999) then
                print*, 'At j = ',j,', i = ',i, 'time =',t, 'interpolation of conveyence was not possible, wl', &
                newY(i,j), 'z',z(i,j),'previous wl',newY(i+1,j), 'previous z',z(i+1,j), 'dimensionless_D(i,j)', &
                dimensionless_D(i,j)
                stop
            end if
			co(i) =q_sk_multi * co(i)
            !if (j .eq. 6) print*, i,'3'


            call r_interpol(elevTable,areaTable,nel,xt,newArea(i,j))
            !if (j .eq. 6) print*, i,'4'
            if (newArea(i,j) .eq. -9999) then
                print*, 'At j = ',j,', i = ',i, 'time =',t, 'interpolation of newArea was not possible'
                stop
            end if
            call r_interpol(elevTable,pereTable,nel,xt,pere(i,j))
            !if (j .eq. 6) print*, i,'5'
            call r_interpol(elevTable,topwTable,nel,xt,bo(i,j))
            call r_interpol(elevTable,skkkTable,nel,xt,sk(i,j))
            !if (j .eq. 6) print*, i,'6'

            !sfi = ( qp(i,j) / co(i) ) ** 2.0
            sfi = qp(i,j) * abs(qp(i,j)) / ( co(i)** 2.0 )

            ! Celerity actual equation:
            !!! new for dkdh
            !do pp = 2,nel
            !    if (newY(i,j) .le. elevTable(pp)) then
            !        dkdh  =(convTable(pp)-convTable(pp-1))/(elevTable(pp)-elevTable(pp-1))
            !        EXIT
            !    endif
            !    if (pp .eq. nel) dkdh =(convTable(pp)-convTable(pp-1))/(elevTable(pp)-elevTable(pp-1))
            !end do
            !if (i .gt. 1) then
            !    dbdx(i)=(bo(i,j)-bo(i-1,j))/dx(i-1,j)
            !else
            !    dbdx(i)=dbdx(i+1)
            !end if

            ! Diffusivity actual equation:
            !diffusivity2(i) = co(i) * co(i) / 2.0 / qp(i,j) / bo(i,j)

            !celerity2(i) = qp(i,j) / co(i) * dkdh / bo(i,j) !+ dbdx(i) * diffusivity2(i) / bo(i,j)


            ! Calculating empirical celerity:
            !celerity2(i)=5.0 / 3.0 * sfi ** 0.3 * abs(qp(i,j)) ** 0.4 / bo(i,j) ** 0.4 / (1./(sk(i,j)*q_sk_multi)) ** 0.6
            !celerity2(i)=5.0 / 3.0 * abs(sfi) ** 0.3 * abs(qp(i,j)) ** 0.4 / bo(i,j) ** 0.4 / (1./(sk(i,j)*q_sk_multi)) ** 0.6

             !min(rightBank(i,j)-leftBank(i,j), bo(i,j))
            chnWidth = rightBank(i,j)-leftBank(i,j)
            chnWidth = min(chnWidth,bo(i,j))
            !celerity2(i)=5.0 / 3.0 * abs(sfi) ** 0.3 * abs(qp(i,j)) ** 0.4 / chnWidth ** 0.4 / (1./(sk(i,j)*q_sk_multi)) ** 0.6

            if (depthCalOk(i) .eq. 1) then
                celerity2(i)=5.0 / 3.0 * abs(sfi) ** 0.3 * abs(qp(i,j)) ** 0.4 / bo(i,j) ** 0.4 / (1./(sk(i,j)*q_sk_multi)) ** 0.6
                diffusivity2(i) = abs(qp(i,j)) / 2.0 / bo(i,j) / abs(sfi)
                vel = qp(i,j)/newArea(i,j)
                if (celerity2(i) .gt. 3.0*vel) celerity2(i) = vel*3.0
            else
                if (qp(i,j) .lt. 1) then
                    celerity2(i)=0.5
                else
                    celerity2(i)=1.0
                end if

                diffusivity2(i)=diffusivity(i,j)
            end if

            !celerity2 = 1.0
            !diffusivity2 = 1000.0



            !if (qp(i,j) .le. 0.05) celerity2(i) = vel       ! change 20210420
            !celerity2(i)=max(celerity2(i),1.25*vel)
            !print*, j, i,newY(i,j),sk(i,j), bo(i,j), celerity2(i)

            !if (j .eq. 6) print*, i, 'sfi=', sfi, 'new depth', newY(i,j) - z(i,j)
            !diffusivity2(i) = abs(qp(i,j)) / 2.0 / bo(i,j) / sfi

            !if (qp(i,j) .le. 0.05) diffusivity2(i) = 1.0       ! change 20210420



        newMassBalance =0
        if (newMassBalance .eq. 1) then
            if (i .gt. 1) then
                !print*, j, i, qp(i-1,j), qp(i,j),newArea(i,j),oldArea(i,j),oldArea(i-1,j)
                !pause 0005
                newArea(i-1,j) = oldArea(i-1,j) + oldArea(i,j) - newArea(i,j) - 2.*dtini/dx(i-1,j)*(qp(i,j)-qp(i-1,j))  ! change 20210407

                !if (j .eq. 6) print*, i,'7',newArea(i-1,j),oldArea(i-1,j),dtini,dx(i-1,j),qp(i,j),qp(i-1,j)
                elevTable = xsec_tab(1,:,i-1,j)
                areaTable = xsec_tab(2,:,i-1,j)
                if ( newArea(i-1,j) .le. 0) then
                    slope = (z(i-1,j)-z(i,j))/dx(i-1,j)
                    call normal_y(i-1, j, q_sk_multi, slope, qp(i-1,j), newY(i-1,j), temp, newArea(i-1,j), temp)
                    currentRoutingNormal(i-1,j) = 1
                    !if ((j .eq. 6) .and. (i .le. 2)) print*, '8 normal'
                else
                    call r_interpol(areaTable,elevTable,nel,newArea(i-1,j),newY(i-1,j))
                    currentRoutingNormal(i-1,j) = 0
                end if
            !if (j .eq. 6) print*, i,'7'
            end if
        else

            if (i .gt. 1) then
           !     print*, 'backward j', j,'i', i, 'routingNotChanged(i-1,j)',routingNotChanged(i-1,j), &
           !     'currentRoutingNormal(i-1,j)',currentRoutingNormal(i-1,j),'dimensionless_D(i-1,j)',dimensionless_D(i-1,j)
                !! If routing method is changed just a few time steps ago, we maintain the same routing to avoid oscillation
                if ( (routingNotChanged(i-1,j) .lt. minNotSwitchRouting2) .and. (currentRoutingNormal(i-1,j) .lt. 3) ) then
                    !print*, routingNotChanged(i-1,j), minNotSwitchRouting2, currentRoutingNormal(i-1,j)
                    !pause 2020
                    if (currentRoutingNormal(i-1,j) .eq. 0) then
                        !newY(i-1,j) = newY(i,j) + sign ( sfi, qp(i,j) ) * dx(i-1,j)
                        newY(i-1,j) = newY(i,j) + sfi * dx(i-1,j)
                    else if (currentRoutingNormal(i-1,j) .eq. 1) then
                        slope = (z(i-1,j)-z(i,j))/dx(i-1,j)
                        if (slope .le. 0.0001) slope = 0.0001
                        call calc_q_sk_multi(i-1,j,qp(i-1,j),q_sk_multi)
                        ! applying normal depth to all the nodes
                        !print*, 'check 1',j,i
                        call normal_crit_y(i-1, j, q_sk_multi, slope, qp(i-1,j), newY(i-1,j), temp, newArea(i-1,j), temp)
                    end if
                else
                    !print*, t/60.-t0, 'hour'
                    !pause 1900
                    !! If DSP: D is below 1.0, we switch to partial diffusive routing
                    if (dimensionless_D(i-1,j) .lt. D_lim1) then
                        slope = (z(i-1,j)-z(i,j))/dx(i-1,j)
                        if (slope .le. 0.0001) slope = 0.0001
                        call calc_q_sk_multi(i-1,j,qp(i-1,j),q_sk_multi)
                        ! applying normal depth to all the nodes
                        call normal_crit_y(i-1, j, q_sk_multi, slope, qp(i-1,j), newY(i-1,j), temp, newArea(i-1,j), temp)
                        !if (slope .eq. 0) slope = TOLERANCE
                         ! Book-keeping: changing from full diffusive to partial diffusive
                        if ( currentRoutingNormal(i-1,j) .ne. 1 ) routingNotChanged(i-1,j) = 0
                        currentRoutingNormal(i-1,j) = 1
                        !print*, 'check 2-1',j,i

                    !pause 1000
                    !! If DSP: D is not below 1.0, we switch to full diffusive routing
                    elseif ( (dimensionless_D(i-1,j) .ge. D_lim1) .and. (dimensionless_D(i-1,j) .lt. D_lim2) ) then
                        print*, 'partial diffusive at j', j, 'i-1',i-1
                        slope = (z(i-1,j)-z(i,j))/dx(i-1,j)
                        if (slope .le. 0.0001) slope = 0.0001
                        call calc_q_sk_multi(i-1,j,qp(i-1,j),q_sk_multi)
                        ! applying normal depth to all the nodes
                        !print*, 'check 3',j,i
                        call normal_crit_y(i-1, j, q_sk_multi, slope, qp(i-1,j), stg1, temp, newArea(i-1,j), temp)

                        !stg2 = newY(i,j) + sign ( sfi, qp(i,j) ) * dx(i-1,j)
                        stg2 = newY(i,j) + sfi * dx(i-1,j)
                        newY(i-1,j) = ( stg2 * (dimensionless_D(i-1,j) - D_lim1) + &
                                        stg1 * (D_lim2 - dimensionless_D(i-1,j)) ) / (D_lim2 - D_lim1)
                         ! Book-keeping: changing from full diffusive to partial diffusive
                        if ( currentRoutingNormal(i-1,j) .ne. 3 ) routingNotChanged(i-1,j) = 0
                        currentRoutingNormal(i-1,j) = 3

                        !pause 1000

                    else
                        !newY(i-1,j) = newY(i,j) + sfi * dx(i-1,j)
                        ! implicit calculation of y
!                        if ( qp(i-1,j) .le. 0.05 ) then
!                            slope = (z(i-1,j)-z(i,j))/dx(i-1,j)
!                            if (slope .le. 0.0001) slope = 0.0001
!                            call calc_q_sk_multi(i-1,j,qp(i-1,j),q_sk_multi)
!                            call normal_crit_y(i-1, j, q_sk_multi, slope, abs(qp(i-1,j)), newY(i-1,j), temp, newArea(i-1,j), temp)
!                        else


                        slope = (z(i-1,j)-z(i,j))/dx(i-1,j)
                        depthYi = newY(i,j) - z(i,j)

                        tempDepthi_1 = oldY(i-1,j)-z(i-1,j)

                        elevTable_1 = xsec_tab(1,:,i-1,j)
                        areaTable_1 = xsec_tab(2,:,i-1,j)
                        pereTable_1 = xsec_tab(3,:,i-1,j)
                        rediTable_1 = xsec_tab(4,:,i-1,j)
                        convTable_1 = xsec_tab(5,:,i-1,j)
                        topwTable_1 = xsec_tab(6,:,i-1,j)
                        skkkTable_1 = xsec_tab(11,:,i-1,j)


                        !tempDepthi_1 = elevTable_1(100) - z(i-1,j)

                        currentSquaredDepth_1=(elevTable_1-z(i-1,j))**2.

                        toll = 1.0
                        iii = 0

!                        ! Applying Newton–Raphson method
                        if (newtonRaphson .eq. 1) then
                        do while ( abs(toll) .gt. 0.001)
                            iii = iii +1
                            tempY_1 = tempDepthi_1 + z(i-1,j)

                            call r_interpol(currentSquaredDepth_1,convTable_1,nel,(tempDepthi_1)**2.0,tempCo_1)
                            call calc_q_sk_multi(i-1,j,qp(i-1,j),temp_q_sk_multi_1)
                            tempCo_1 = tempCo_1 * temp_q_sk_multi_1

                            call r_interpol(elevTable_1,areaTable_1,nel,tempY_1,tempArea_1)
                            call r_interpol(elevTable_1,pereTable_1,nel,tempY_1,tempPere_1)
                            call r_interpol(elevTable_1,rediTable_1,nel,tempY_1,tempRadi_1)
                            call r_interpol(elevTable_1,topwTable_1,nel,tempY_1,tempbo_1)
                            call r_interpol(elevTable_1,skkkTable_1,nel,tempY_1,tempsk_1)!

                            tempsfi_1 = qp(i-1,j) * abs(qp(i-1,j)) / ( tempCo_1** 2.0 )!

                            ffy = tempDepthi_1 - depthYi + dx(i-1,j)*slope - 0.5*dx(i-1,j)*(sfi + tempsfi_1)

                            dkda=tempsk_1*((5.0/3.0*tempArea_1**(2.0/3.0)*tempPere_1)- &
                                (tempArea_1**(5.0/3.0)*2.0/tempbo_1))/(tempPere_1**(5.0/3.0))

                            ffprime = 1 + dx(i-1,j) * tempbo_1 *  qp(i-1,j) * abs(qp(i-1,j)) / (tempCo_1 ** 3.0) * dkda

                            tempDepthi_1_new = tempDepthi_1 - ffy / ffprime

                            tempDepthi_1_new = max(tempDepthi_1_new,0.005)

                            toll = abs(tempDepthi_1_new - tempDepthi_1)

                            !if ( (j .eq. 3) .and. (i-1 .eq. 3) ) then
                            !    if ((t+dtini/60.-1267)*(t+dtini/60.-1275) .le. 0) then
                            !        print*, t+dtini/60., j, i-1, iii, qp(i-1,j), tempCo_1, sfi, tempsfi_1, ffy, ffprime, &
                            !        tempDepthi_1, tempDepthi_1_new
                            !        !pause
                            !    end if
                            !end if

                            if(iii .gt. 20)then
                                print*, 'Warning: Depth iteration reached maximum trial at j=', j, 'i=', i-1 , 'and',i, &
                                'depths are',tempDepthi_1, tempDepthi_1_new, 'slope=', slope, 'tempbo_1=',tempbo_1, 'dx=', dx(i-1,j)
                                print*, 'depth at d/s',depthYi
                                !tempDepthi_1 = oldY(i-1,j)-z(i-1,j)
                                depthCalOk(i-1) = 0
                                EXIT
                            endif
                            tempDepthi_1 = tempDepthi_1_new
                            depthCalOk(i-1) = 1



                        end do
                        end if

                        ! Applying mid point bisection
                        if (newtonRaphson .eq. 0) then
                        tempY_1 = elevTable_1(2)
                        tempY_2 = depthYi * 3. + z(i-1,j)
                        tempY_3 = (tempY_1 + tempY_2) / 2.
                        do while ( abs(toll) .gt. 0.001)
                            iii = iii +1

                            call r_interpol(currentSquaredDepth_1,convTable_1,nel,(tempY_1-z(i-1,j))**2.0,tempCo_1)
                            call r_interpol(currentSquaredDepth_1,convTable_1,nel,(tempY_2-z(i-1,j))**2.0,tempCo_2)
                            call r_interpol(currentSquaredDepth_1,convTable_1,nel,(tempY_3-z(i-1,j))**2.0,tempCo_3)

                            call calc_q_sk_multi(i-1,j,qp(i-1,j),temp_q_sk_multi_1)
                            tempCo_1 = tempCo_1 * temp_q_sk_multi_1
                            tempCo_2 = tempCo_2 * temp_q_sk_multi_1
                            tempCo_3 = tempCo_3 * temp_q_sk_multi_1

                            tempsfi_1 = qp(i-1,j) * abs(qp(i-1,j)) / ( tempCo_1** 2.0 )
                            tempsfi_2 = qp(i-1,j) * abs(qp(i-1,j)) / ( tempCo_2** 2.0 )
                            tempsfi_3 = qp(i-1,j) * abs(qp(i-1,j)) / ( tempCo_3** 2.0 )

                            ffy_1 = (tempY_1-z(i-1,j)) - depthYi + dx(i-1,j)*slope - 0.5*dx(i-1,j)*(sfi + tempsfi_1)
                            ffy_2 = (tempY_2-z(i-1,j)) - depthYi + dx(i-1,j)*slope - 0.5*dx(i-1,j)*(sfi + tempsfi_2)
                            ffy_3 = (tempY_3-z(i-1,j)) - depthYi + dx(i-1,j)*slope - 0.5*dx(i-1,j)*(sfi + tempsfi_3)

                            if ((ffy_1 * ffy_2) .gt. 0.) then
                                tempY_2 = (tempY_2 - z(i-1,j)) * 2.0 + z(i-1,j)
                            elseif ((ffy_1 * ffy_3) .le. 0.) then
                                tempY_2 = tempY_3
                            elseif ((ffy_2 * ffy_3) .le. 0.) then
                                tempY_1 = tempY_3
                            end if
                            tempY_3 = (tempY_1 + tempY_2) / 2.0
                            toll = tempY_2 - tempY_1
                            tempDepthi_1 = tempY_3 - z(i-1,j)
                            depthCalOk(i-1) = 1


                        end do
                        end if      ! calculating depth based on d/s WL starts

                        newY(i-1,j) = tempDepthi_1 + z(i-1,j)

                        call r_interpol(elevTable_1,areaTable_1,nel,tempY_1,tempArea_1)
                        call r_interpol(elevTable_1,topwTable_1,nel,tempY_1,tempbo_1)

                        usFroud = abs(qp(i-1,j))/sqrt(grav*tempArea_1**3.0/tempbo_1)

                        !if (usFroud .ge. 1.0) then
                        !   call normal_crit_y(i-1, j, q_sk_multi, slope, qp(i-1,j), newY(i-1,j), temp, newArea(i-1,j), temp)
                        !   depthCalOk(i-1) = 0
                        !end if

                        if (newY(i-1,j) .gt. 10.0**5.) newY(i-1,j) = 10.0**5.
!                        end if





                         ! Book-keeping: changing from partial diffusive to full diffusive
                        if ( currentRoutingNormal(i-1,j) .ne. 0 ) routingNotChanged(i-1,j) = 0
                        currentRoutingNormal(i-1,j) = 0

                    end if
                end if

                if (newY(i-1,j)-z(i-1,j) .le. 0.) then
                    print*, 'depth is negative at time=,', t,'j= ', j,'i=',i-1,'newY=',(newY(jj,j),jj=1,ncomp)
                    print*, 'dimensionless_D',(dimensionless_D(jj,j),jj=1,ncomp)
                    print*, 'newQ',(newQ(jj,j),jj=1,ncomp)
                    print*, 'Bed',(z(jj,j),jj=1,ncomp)
                    print*, 'dx',(dx(jj,j),jj=1,ncomp-1)
                    pause 777
                end if

            end if
        end if

                ! Book-keeping: Counting the number as for how many time steps the routing method is unchanged
                routingNotChanged(i-1,j) = routingNotChanged(i-1,j) + 1

        end do

            !celerity(i) = sign ( sum(celerity2) / ncomp, qp(i) )


           ! if (sum(depthCalOk) .lt. ncomp) then
           !     print*, 'celerity 2', celerity2(1:ncomp)
           !     print*, 'celerity 1', celerity(1:ncomp,j)
           ! end if



            celerity(1:ncomp,j) =  sum(celerity2(1:ncomp)) / ncomp ! change 20210524
            !celerity(1:ncomp,j) = minval(celerity2(1:ncomp)) ! change 20210423


            if (celerity(1,j) .lt. 0.5) celerity(1:ncomp,j) = 0.5
            !celerity(1:ncomp,j) =  celerity2(1:ncomp)


        !print*,j, diffusivity2(1:ncomp)
        diffusivity(1:ncomp,j)=sum(diffusivity2(1:ncomp)) / ncomp
        !diffusivity = 1000.
		do i = 1, ncomp
			if (diffusivity(i,j) .gt. maxDiffuLm) diffusivity(i,j) = maxDiffuLm !!! Test
			if (diffusivity(i,j) .lt. minDiffuLm) diffusivity(i,j) = minDiffuLm !!! Test
		end do

		!diffusivity = 200.         ! Test
		!celerity = 1.              ! Test



!!++++++++++++++++++++ Diffusive wave Backward sweep ends +++++++++++++++++++!!




end subroutine mesh_diffusive_backward
