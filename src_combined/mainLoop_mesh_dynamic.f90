subroutine mesh_dynamic(dtini_given, ppp,qqq, t0, t, tfin, saveInterval)

    use constants_module
    use arrays_module
    use matrix_module
    use var_module
    use arrays_section_module
    use xsec_attribute_module
    use sgate_module

    implicit none

    integer(kind=4), intent(in) :: ppp, qqq
    real(kind=4), intent(in) :: dtini_given, t0, t, tfin, saveInterval

    integer(kind=4) :: i, j, n, ntim, igate, pp, boundaryFileMaxEntry, noLatFlow, noQSKtable, saveFrequency
    real(kind=4) :: cour, da, dq, dxini, yy, x, thetas, thesinv
    real(kind=4) :: skk, qq, qn, xt, r_interpol, maxCourant, sfi, r_interpo_nn, frds
    real(kind=4) :: arean, areac, hyrdn, hyrdc, perimn, perimc, qcrit, s0ds, timesDepth, latFlowValue, q_sk_multi



    dtini = lastKnownDiffuDT

    !print*, 'in dynamic', dtini, dtini_given

		call correctDT(t0, t, saveInterval, tfin)
    !print*, 'Dynamic dtini', dtini
    !pause 5000
            ! Set upstream discharge
        dqp(1) = newQ(1) - oldQ(1)
        dap(1) = 0.0
        !print*, 'areas', area_norm, oldArea(1), newArea(1)


        !temp disabled!
        if (S_0 .gt. 0) then
            !print*, y_crit_us-z(1), y_norm_us-z(1)!; pause 5000
            if(y_crit_us .ge. y_norm_us) then ! Fr_us>1: supercritical flow
              if (frus2 .ge. 1.0) then
                print*, 'US boundary is supercritical'
                newY(1) = y_norm_us
                newArea(1) = area_norm
                dap(1)=area_norm -  oldArea(1)   !update from boundary condition time series
              endif
            endif
        endif
        !print*, y_norm_us, y_crit_us, dap(1)
        !pause 5000



        call section()
        ! Nazmul: The subroutine calls the attribute tables and interpolate according to the available water level
        thes=thetas

        call matrixp()

        do i=2,ncomp
            !cour=dt(i)/dx(i-1)
            cour=dtini/dx(i-1)
            !rhs1=-cour*(f1(i)-f1(i-1)-d1(i)+d1(i-1))+lateralFlow(i-1)*dtini !*dx(i)/dx(i-1)
            rhs1=-cour*(f1(i)-f1(i-1)-d1(i)+d1(i-1))+lateralFlow(i)*dtini !*dx(i)/dx(i-1)
            rhs2=-cour*(f2(i)-f2(i-1)-d2(i)+d2(i-1))+dtini*grav*(ci2(i)+aso(i))
            c11=g11inv(i)*b11(i-1)+g12inv(i)*b21(i-1)
            c12=g11inv(i)*b12(i-1)+g12inv(i)*b22(i-1)
            c21=g21inv(i)*b11(i-1)+g22inv(i)*b21(i-1)
            c22=g21inv(i)*b12(i-1)+g22inv(i)*b22(i-1)
            dap(i)=g11inv(i)*rhs1+g12inv(i)*rhs2-c11*dap(i-1)-c12*dqp(i-1)
            dqp(i)=g21inv(i)*rhs1+g22inv(i)*rhs2-c21*dap(i-1)-c22*dqp(i-1)
           ! print*,'rhs1=', rhs1, rhs2, c11, c12, c21, c22
        end do


!        if (option_dsbc.eq.1) then
!            dac(ncomp)=dap(ncomp)
!            yn=(area(ncomp)+dap(ncomp))/bo(ncomp)
!            arean=yn*bo(ncomp)
!            perimn=2.0*yn+bo(ncomp)
!            hyrdn=arean/perimn
!            s0ds=-((z(ncomp)-z(ncomp-1))/dx(ncomp))
!            qn=skk*arean*hyrdn**(2.0/3.0)*sqrt(s0ds)
!            dqp(ncomp)=qn-oldQ(ncomp)
!            dqc(ncomp)=dqp(ncomp)
!
!       elseif(option_dsbc.eq.2)then
!            dac(ncomp)=dap(ncomp)
!            yn=(area(ncomp)+dap(ncomp))/bo(ncomp)
!            areac=yn*bo(ncomp)
!            perimc=2.0*yn+bo(ncomp)
!            hyrdc=areac/perimc
!            s0ds=-((z(ncomp)-z(ncomp-1))/dx(ncomp))
!            qn=skk*areac*hyrdc**(2.0/3.0)*sqrt(s0ds)
!            qcrit=1.05*(((yn**3.0)*(bo(ncomp)**2.0)*grav)**(1.0/2.0))
!            write(*,*)qcrit
!            dqp(ncomp)=qcrit-oldQ(ncomp)
!            dqc(ncomp)=dqp(ncomp)

 !       else
            !dac(ncomp)=0.0
            !dap(ncomp)=0.0
            !dqc(ncomp)=dqp(ncomp)

!            dap(ncomp)=0.0	!checked email !for critical comment out
! change for unsteady flow
			!dap(ncomp) = newArea(ncomp) - oldArea(ncomp)
            !dac(ncomp)=dap(ncomp)	!checked email
            !dqc(ncomp)=dqp(ncomp)	!checked email
            !dac(ncomp)=newArea(ncomp) - oldArea(ncomp)
 !       endif

        !! New addition 20200624 start
        !Check Froud number at the downstream end


        frds=abs(oldQ(ncomp) + dqp(ncomp))/sqrt(grav*(newArea(ncomp))**3.0/bo(ncomp))   !!! applied q(n,ncomp)

       !frds=abs(q(n,ncomp) + dqp(ncomp))/sqrt(grav*ads**3.0/bo(ncomp))   !!! applied q(n,ncomp)

       if(frds .lt. 1.0) then
! Fr_ds <1: subcritical flow
!
          dap(ncomp)=newArea(ncomp) - oldArea(ncomp)    !update from downstream time series
          dac(ncomp)=dap(ncomp)
          dqc(ncomp)=dqp(ncomp)
       else
!
! Fr_ds>=1: supercritical
        print*, 'DS boundary is supercritical'
        print*, frds, newY(ncomp)-z(ncomp), oldQ(ncomp), dqp(ncomp), newArea(ncomp), bo(ncomp)!; pause 5000
         dac(ncomp)=dap(ncomp)
         dqc(ncomp)=dqp(ncomp)
       endif

!! New addition 20200624 end

        ! Update via predictor
        areap = area + dap
        qp = oldQ + dqp

        !do i=1,ncomp
        !    if (areap(i) .le. 0.) areap(i) = 0.01
        !end do


        ! applying lateral flow at the qp
        !do i=1,noLatFlow
        !    latFlowValue = r_interpol(lateralFlowTable(1, :, i),lateralFlowTable(2, :, i),dataInEachLatFlow(i),t)
        !    qp(latFlowLocations(i)) = qp(latFlowLocations(i)) + &
        !        latFlowValue*(dx(latFlowLocations(i)-1)+dx(latFlowLocations(i)))*0.5
        !    print*, 'lat flow=', latFlowValue, latFlowLocations(i), &
        !        latFlowValue*(dx(latFlowLocations(i)-1)+dx(latFlowLocations(i)))*0.5, &
        !        qp(latFlowLocations(i))
        !end do

        !print*, 'areap', (areap(i), i=1, ncomp)


        call secpred()
        thes=thesinv
        call matrixc()

        do i=ncomp-1,1,-1
            !cour=dt(i)/dx(i)
            cour=dtini/dx(i)
            !rhs1=-cour*(f1(i+1)-f1(i)-d1(i+1)+d1(i))+lateralFlow(i)*dtini
            rhs1=-cour*(f1(i+1)-f1(i)-d1(i+1)+d1(i))+lateralFlow(i+1)*dtini
            rhs2=-cour*(f2(i+1)-f2(i)-d2(i+1)+d2(i))+dt(i)*grav*(ci2(i)+aso(i))
            c11=g11inv(i)*b11(i+1)+g12inv(i)*b21(i+1)
            c12=g11inv(i)*b12(i+1)+g12inv(i)*b22(i+1)
            c21=g21inv(i)*b11(i+1)+g22inv(i)*b21(i+1)
            c22=g21inv(i)*b12(i+1)+g22inv(i)*b22(i+1)
            dac(i)=g11inv(i)*rhs1+g12inv(i)*rhs2-c11*dac(i+1)-c12*dqc(i+1)
            dqc(i)=g21inv(i)*rhs1+g22inv(i)*rhs2-c21*dac(i+1)-c22*dqc(i+1)
        end do
        ! Upstream boundary condition
        ! Prescribed discharge at the upstream
        ! Area correction is calculated



        if ( S_0 .gt. 0.) then
          if( y_crit_us .ge. y_norm_us ) then ! Steep slope
            frus2 = abs(newQ(1))/sqrt(grav*( dac(1)+areap(1) )**3.0/bo(1))
            !print*, 'frus2', frus2, newQ(1), oldArea(1), areap(1), dac(1), bo(1)
            !pause 5000
            if (frus2 .ge. 1) then   ! Fr_us>=1: supercritical flow
              dqc(1)=dqp(1)
              dac(1)=dap(1)
            else
              dqc(1)=dqp(1)	!checked email
              dap(1)=dac(1)	!checked email
            endif
          else
            dqc(1)=dqp(1)	!checked email
            dap(1)=dac(1)	!checked email
          endif
        else
          dqc(1)=dqp(1)	!checked email
          dap(1)=dac(1)	!checked email
        endif



        ! Final update
        do i=1,ncomp
            da=(dap(i)+dac(i))/2.0
            dq=(dqp(i)+dqc(i))/2.0
            newArea(i)=da+area(i)
            if(newArea(i) <= 0.0) newArea(i)=0.001

!           Now calculate y based on area calculated
!-------------------------------------
            elevTable(:) = xsec_tab(1,:,i)
            areaTable(:) = xsec_tab(2,:,i)

    !       interpolate the cross section attributes based on FINAL CALCULATED area
            xt=newArea(i)
            newY(i)=r_interpol(areaTable,elevTable,nel,xt)
!-------------------------------------

            newQ(i)=oldQ(i)+dq
            !froud(i)=abs(newQ(i))/sqrt(grav*newArea(i)**3.0/bo(i))

        end do
        !print*, newArea(1)
        !pause 5001

        do i=1,ncomp
            q_sk_multi = 1.0
            do pp = 1, size(Q_sk_tableEntry)
                if (  ( eachQSKtableNodeRange(1,pp) - i) * ( eachQSKtableNodeRange(2,pp) - i) .le. 0 ) then
                    q_sk_multi = r_interpo_nn(Q_Sk_Table(1,1:Q_sk_tableEntry(pp),pp),   &
                      Q_Sk_Table(2,1:Q_sk_tableEntry(pp),pp),Q_sk_tableEntry(pp),qp(i))
                end if
            end do

            sfi = ( newQ(i) / co(i) ) ** 2.0
            celerity2(i)=5.0 / 3.0 * sfi ** 0.3 * abs(newQ(i)) ** 0.4 / bo(i) ** 0.4 / (1/(sk(i)*q_sk_multi)) ** 0.6
        end do
        celerity =  sum(celerity2) / ncomp



end subroutine mesh_dynamic
