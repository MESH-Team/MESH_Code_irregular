subroutine mesh_dynamic_predictor(dtini_given, t0, t, tfin, saveInterval,j)

    use constants_module
    use arrays_module
    use matrix_module
    use var_module
    use arrays_section_module
    use xsec_attribute_module
    use sgate_module

    implicit none

    integer(kind=4), intent(in) :: j
    real(kind=4), intent(in) :: dtini_given, t0, t, tfin, saveInterval

    integer(kind=4) :: i, igate, pp, k, nodenb, linknb
    real(kind=4) :: cour, frds, areasum, yk_ncomp, yav, areak_ncomp, areav
    real(kind=4) :: xt, r_interpol, r_interpo_nn
    real(kind=4) :: q_sk_multi



    !dtini = dtini_given

    !print*, 'in dynamic', dtini, dtini_given


    !cannot calculate dt inside this subroutine
		!call correctDT(t0, t, saveInterval, tfin)

    !print*, 'Dynamic dtini', dtini
    !pause 5000
        dtini = dtini_given

        !print*, 'done'; pause 5001

        call section(j)



        !print*, 'areas', area_norm, oldArea(1), newArea(1)
        if (S_0 .gt. 0) then
            if(y_crit_us .ge. y_norm_us) then ! Fr_us>1: supercritical flow
              if (frus2 .ge. 1.0) then
                print*, 'US boundary is supercritical '!; pause 5000
                newY(1,j) = y_norm_us
                newArea(1,j) = area_norm
                dap(1,j)=area_norm -  oldArea(1,j)   !update from boundary condition time series
              endif
            endif
        endif
        !print*, y_norm_us, y_crit_us, dap(1)
        !pause 5000







            !+++----------------------------------------------------------------
            !+ Hand over water from upstream to downstream properly according
            !+ to the nature of link connections, i.e., serial or branching.
            !+ Refer to p.52,RM1_MESH
            !+++----------------------------------------------------------------
            if (ndep(j).gt.0) then
                !*total water areas at n+1 at the end nodes of upstream links that join link j
                areasum=0.0
                do k=1, ndep(j)
                    linknb=uslinks(k,j); nodenb=nx1(linknb)
                    areasum=areasum + oldArea(nodenb,linknb) + dap(nodenb,linknb)
                end do

                dqp(1,j)=0.0;

    print*, 'Nazmul'
                yav=0.0
                do k=1, ndep(j)
                    linknb=uslinks(k,j); nodenb=nx1(linknb)
                    !**dqp(1,j)
                    dqp(1,j)=dqp(1,j)+dqp(nodenb,linknb)
                    !**dap(1,j)
                    !*area at the end nod of link k at time n+1
                    areak_ncomp = oldArea(nodenb,linknb) + dap(nodenb,linknb)
                    !bok_ncomp=bo(nodenb,linknb) !* bottom width at the end of link k in the immed. upstream of link j.
                    !sslp=traps(nodenb, linknb)  !* side slope of trapezoidal channel.

                    elevTable = xsec_tab(1,:,nodenb,linknb)
                    areaTable = xsec_tab(2,:,nodenb,linknb)
                    yk_ncomp = r_interpol(areaTable,elevTable,nel,areak_ncomp)
                    !* weighted average based on areas at the end nodes of upstream link ks
                    yav = yav + (areak_ncomp/areasum)*yk_ncomp  !! Nazmul : Learn this
                end do
                !* Area estimated for time n+1
                elevTable = xsec_tab(1,:,1,j)
                areaTable = xsec_tab(2,:,1,j)
                areav = r_interpol(elevTable,areaTable,nel,yav)
                dap(1,j) = areav - oldArea(1,j)
            else
                ! Set upstream discharge
                dqp(1,j) = newQ(1,j) - oldQ(1,j)
                dap(1,j) = 0.0
            end if

    print*, 'Nazmul'






















        ! Nazmul: The subroutine calls the attribute tables and interpolate according to the available water level
        thes=thetas

        call matrixp(j)

        do i=2,ncomp
            !cour=dt(i)/dx(i-1)
            cour=dtini/dx(i-1,j)
            !rhs1=-cour*(f1(i)-f1(i-1)-d1(i)+d1(i-1))+lateralFlow(i-1)*dtini !*dx(i)/dx(i-1)
            rhs1=-cour*(f1(i)-f1(i-1)-d1(i)+d1(i-1))+lateralFlow(i)*dtini !*dx(i)/dx(i-1)
            rhs2=-cour*(f2(i)-f2(i-1)-d2(i)+d2(i-1))+dtini*grav*(ci2(i)+aso(i))
            c11=g11inv(i)*b11(i-1)+g12inv(i)*b21(i-1)
            c12=g11inv(i)*b12(i-1)+g12inv(i)*b22(i-1)
            c21=g21inv(i)*b11(i-1)+g22inv(i)*b21(i-1)
            c22=g21inv(i)*b12(i-1)+g22inv(i)*b22(i-1)
            dap(i,j)=g11inv(i)*rhs1+g12inv(i)*rhs2-c11*dap(i-1,j)-c12*dqp(i-1,j)
            dqp(i,j)=g21inv(i)*rhs1+g22inv(i)*rhs2-c21*dap(i-1,j)-c22*dqp(i-1,j)
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


        frds=abs(oldQ(ncomp,j) + dqp(ncomp,j))/sqrt(grav*(newArea(ncomp,j))**3.0/bo(ncomp))   !!! applied q(n,ncomp)

       !frds=abs(q(n,ncomp) + dqp(ncomp))/sqrt(grav*ads**3.0/bo(ncomp))   !!! applied q(n,ncomp)

       if(frds .lt. 1.0) then
! Fr_ds <1: subcritical flow
!
          dap(ncomp,j)=newArea(ncomp,j) - oldArea(ncomp,j)    !update from downstream time series
          dac(ncomp,j)=dap(ncomp,j)
          dqc(ncomp,j)=dqp(ncomp,j)
       else


!
! Fr_ds>=1: supercritical
        print*, 'DS boundary is supercritical'
        !print*, frds, y(n+1,ncomp)!; pause 5000
         dac(ncomp,j)=dap(ncomp,j)
         dqc(ncomp,j)=dqp(ncomp,j)
       endif


!! New addition 20200624 end

        ! Update via predictor
        areap(:,j) = area(:) + dap(:,j)
        qp(:,j) = oldQ(:,j) + dqp(:,j)

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


        !+++-------------CHECK: HOW IT WORKS-------------------------------------------------+
        !+ Set downstream boundary conditions of the most downstream link
        !+ eventually for corrector step. It is assumed that stage
        !+ data is always available.
        !+++--------------------------------------------------------------+
        !* Rec. channel: estimate y at the final node and time n
     !   if (j==nlinks) then
     !       call dsbc(n, j)
     !   end if

end subroutine mesh_dynamic_predictor
