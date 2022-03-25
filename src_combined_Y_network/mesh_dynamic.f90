subroutine mesh_dynamic_predictor(dtini_given, t0, t, tfin, saveInterval,j)

    use constants_module
    use arrays_module
    use matrix_module
    use var_module
    use arrays_section_module
    use xsec_attribute_module
    use sgate_module
    use subtools

    implicit none

    integer, intent(in) :: j
    real, intent(in) :: dtini_given, t0, t, tfin, saveInterval

    integer :: i, igate, pp, k, nodenb, linknb
    real :: cour, frds, areasum, yk_ncomp, yav, areak_ncomp, areav, sumOldQ
    real :: xt, r_interpo_nn
    real :: q_sk_multi




        !print*, 'areas', area_norm, y_norm_us, oldArea(1,j), newArea(1,j)
 !       if (S_0 .gt. 0) then
 !           if(y_crit_us .ge. y_norm_us) then ! Fr_us>1: supercritical flow
 !             if (frus2 .ge. 1.0) then
 !               print*, 'US boundary is supercritical at river reach ', j!; pause 5000
 !               newY(1,j) = y_norm_us
 !               newArea(1,j) = area_norm
 !               dap(1,j)=area_norm -  oldArea(1,j)   !update from boundary condition time series
 !             endif
 !           endif
 !       endif
        !print*, y_norm_us, y_crit_us,

!print*, dap(1:ncomp,j); pause
!print*, dqp(1:ncomp,j); pause
!print*, j, 'newQ', newQ(1,j), newArea(ncomp,j)
!print*, 'dqp', dqp(1:ncomp,j), dap(1:ncomp,j)
!pause
        call section(j)

        !print*, j, 'section complete'


        ! Nazmul: The subroutine calls the attribute tables and interpolate according to the available water level
        thes=thetas

        call matrixp(j)

        !print*, j, 'matrixp complete'

        do i=2,ncomp
            !cour=dt(i)/dx(i-1)
            cour=dtini/dx(i-1,j)

            rhs1=-cour*(f1(i)-f1(i-1)-d1(i)+d1(i-1))+lateralFlow(i,j)*dtini     !
            !rhs1=-cour*(f1(i)-f1(i-1)-d1(i)+d1(i-1))+0.5*(lateralFlow(i,j)+lateralFlow(i-1,j)/dx(i-1,j)*dx(i-2,j) )*dtini       ! change 20210315  ! brings difficulties ! reverting
            rhs2=-cour*(f2(i)-f2(i-1)-d2(i)+d2(i-1))+dtini*grav*(ci2(i)+aso(i,j))
            c11=g11inv(i)*b11(i-1)+g12inv(i)*b21(i-1)
            c12=g11inv(i)*b12(i-1)+g12inv(i)*b22(i-1)
            c21=g21inv(i)*b11(i-1)+g22inv(i)*b21(i-1)
            c22=g21inv(i)*b12(i-1)+g22inv(i)*b22(i-1)
            dap(i,j)=g11inv(i)*rhs1+g12inv(i)*rhs2-c11*dap(i-1,j)-c12*dqp(i-1,j)
            dqp(i,j)=g21inv(i)*rhs1+g22inv(i)*rhs2-c21*dap(i-1,j)-c22*dqp(i-1,j)
           ! print*,'rhs1=', rhs1, rhs2, c11, c12, c21, c22
        end do

        !print*, j, dap(:,j)
        !print*, j, dqp(:,j)


        ! Update via predictor
        ! The areap(ncomp,j) is later corrected before the corrector loop starts
        areap(1:ncomp,j) = area(1:ncomp) + dap(1:ncomp,j)
        qp(1:ncomp,j) = oldQ(1:ncomp,j) + dqp(1:ncomp,j)

        !print*, j, qp(1:ncomp,j)
        !print*, 'end of forward', j, lateralFlow(1:ncomp,j)


end subroutine mesh_dynamic_predictor

subroutine mesh_dynamic_corrector(dtini_given, t0, t, tfin, saveInterval,j)

    use constants_module
    use arrays_module
    use matrix_module
    use var_module
    use arrays_section_module
    use xsec_attribute_module
    use sgate_module
    use subtools

    implicit none

    integer, intent(in) :: j
    real, intent(in) :: dtini_given, t0, t, tfin, saveInterval

    integer :: i, k, igate, pp, tableLength, linknb_ds, linknb_us, nodenb
    real :: da, dq, sfi, cour, currentQ,linknb, qnp1_us, qnp1_ds, qsum
    real :: xt, r_interpo_nn
    real :: q_sk_multi


        call secpred(j)

        !print*, j, 'secpred compltete', ncomp

        thes=thesinv
        call matrixc(j)

        !print*, j, 'matrixc compltete', ncomp


        do i=ncomp-1,1,-1
            !cour=dt(i)/dx(i)
            cour=dtini/dx(i,j)
            !rhs1=-cour*(f1(i+1)-f1(i)-d1(i+1)+d1(i))+lateralFlow(i,j)*dtini
            rhs1=-cour*(f1(i+1)-f1(i)-d1(i+1)+d1(i))+lateralFlow(i+1,j)*dtini

            !rhs1=-cour*(f1(i+1)-f1(i)-d1(i+1)+d1(i))+0.5*(lateralFlow(i+1,j)+lateralFlow(i,j)/dx(i,j)*dx(i-1,j))*dtini      ! change 20210315  ! change 20210315  ! brings difficulties ! reverting

            rhs2=-cour*(f2(i+1)-f2(i)-d2(i+1)+d2(i))+dt(i)*grav*(ci2(i)+aso(i,j))
            c11=g11inv(i)*b11(i+1)+g12inv(i)*b21(i+1)
            c12=g11inv(i)*b12(i+1)+g12inv(i)*b22(i+1)
            c21=g21inv(i)*b11(i+1)+g22inv(i)*b21(i+1)
            c22=g21inv(i)*b12(i+1)+g22inv(i)*b22(i+1)
            dac(i,j)=g11inv(i)*rhs1+g12inv(i)*rhs2-c11*dac(i+1,j)-c12*dqc(i+1,j)
            dqc(i,j)=g21inv(i)*rhs1+g22inv(i)*rhs2-c21*dac(i+1,j)-c22*dqc(i+1,j)
        end do
        ! Upstream boundary condition
        ! Prescribed discharge at the upstream
        ! Area correction is calculated
!print*, j, lateralFlow(:,j)
!pause
! ignoring the upstream critical condition check for the large rivers
!        if ( S_0 .gt. 0.) then
!          if( y_crit_us .ge. y_norm_us ) then ! Steep slope
!            frus2 = abs(newQ(1,j))/sqrt(grav*( dac(1,j)+areap(1,j) )**3.0/bo(1,j))
!            print*, 'frus2', frus2, newQ(1,j), oldArea(1,j), areap(1,j), dac(1,j), bo(1,j)
!            !pause 5000
!            if (frus2 .ge. 1) then   ! Fr_us>=1: supercritical flow
!              dqc(1,j)=dqp(1,j)
!              dac(1,j)=dap(1,j)
!            else
!              dqc(1,j)=dqp(1,j)	!checked email
!              dap(1,j)=dac(1,j)	!checked email
!            endif
!          else
!            dqc(1,j)=dqp(1,j)	!checked email
!            dap(1,j)=dac(1,j)	!checked email
!          endif
!        else
!          dqc(1,j)=dqp(1,j)	!checked email
!          dap(1,j)=dac(1,j)	!checked email
!        endif


          dqc(1,j)=dqp(1,j)	!checked email
          dap(1,j)=dac(1,j)	!checked email

        !print*, newArea(35,1),area(35),areap(35,1),'65'

        ! Final update
        do i=1,ncomp

            !print*, j, i, oldArea(i,j), oldY(i,j)-z(i,j)


            da=(dap(i,j)+dac(i,j))/2.0
            dq=(dqp(i,j)+dqc(i,j))/2.0
            newArea(i,j)=da+oldArea(i,j)
            if(newArea(i,j) <= 0.0) newArea(i,j)=0.01

!           Now calculate y based on area calculated
!-------------------------------------
    !       interpolate the cross section attributes based on FINAL CALCULATED area
            xt=newArea(i,j)
            !print*, i, j, ncomp, newArea(35,1),'651'
            if (applyNaturalSection .eq. 0) then
                newY(i,j) = newArea(i,j) / bo(i,j) + z(i,j)
            else
                elevTable = xsec_tab(1,:,i,j)
                areaTable = xsec_tab(2,:,i,j)
                call r_interpol(areaTable,elevTable,nel,xt,newY(i,j))
                if (newY(i,j) .eq. -9999) then
                    print*, 'At j = ',j,', i = ',i, 'time =',t, 'interpolation of newY was not possible'
                    stop
                end if
            end if
            !print*, i, j, ncomp, newArea(35,1),'652'
!-------------------------------------

            newQ(i,j)=oldQ(i,j)+dq
            !froud(i)=abs(newQ(i))/sqrt(grav*newArea(i)**3.0/bo(i))

        end do
        !print*, newArea(1)
        !pause 5001
        !print*, j, dac(1:nx1(j),j)
        !print*, j, dqc(1:nx1(j),j)

        !print*, newArea(35,1),'66'
 !       do i=1,ncomp
 !           currentQ = qp(i,j)
 !           q_sk_multi = 1.0
 !           do pp = 1, noQSKtable(j)
 !               if (  ( eachQSKtableNodeRange(1,pp,j) - i) * ( eachQSKtableNodeRange(2,pp,j) - i) .le. 0 ) then
 !                   tableLength = Q_sk_tableEntry(pp,j)
 !                   q_sk_multi = r_interpo_nn(Q_Sk_Table(1,1:tableLength,pp,j), &
 !                       Q_Sk_Table(2,1:tableLength,pp,j),tableLength,currentQ)
 !               end if
 !           end do
!            co(i) = q_sk_multi*co(i)

 !           sfi = ( newQ(i,j) / co(i) ) ** 2.0
 !           celerity2(i)=5.0 / 3.0 * sfi ** 0.3 * abs(newQ(i,j)) ** 0.4 / bo(i,j) ** 0.4 / (1/(sk(i,j)*q_sk_multi)) ** 0.6
 !       end do
 !       celerity(1:ncomp,j) =  sum(celerity2) / ncomp

end subroutine mesh_dynamic_corrector