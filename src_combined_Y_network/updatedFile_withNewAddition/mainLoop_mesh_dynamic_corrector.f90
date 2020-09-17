subroutine mesh_dynamic_corrector(dtini_given, t0, t, tfin, saveInterval,j)

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

    integer(kind=4) :: i, k, igate, pp, tableLength, linknb_ds, linknb_us, nodenb
    real(kind=4) :: da, dq, sfi, cour, currentQ,linknb, qnp1_us, qnp1_ds, qsum
    real(kind=4) :: xt, r_interpol, r_interpo_nn
    real(kind=4) :: q_sk_multi





        call secpred(j)






            !+++------------------------------------------------------------+
            !+ Handle downstream boundary condition for a link j that has a link
            !+ immediately downstream.
            !+ **Note. dac and dqc at the last node of the most downstream
            !+ link are computed in sub dsbc during predictor step.
            !+ Refer to p.53-1,RM1_MESH
            !+++------------------------------------------------------------+
            if (j<nlinks) then
                linknb=dslink(j)
                !*dac(ncomp,j)

                newArea(1,linknb)=oldArea(1,linknb)+0.5*(dap(1,linknb)+dac(1,linknb))
                xt = newArea(1,linknb)
                elevTable = xsec_tab(1,:,1,linknb)
                areaTable = xsec_tab(2,:,1,linknb)
                newY(1,linknb)=r_interpol(areaTable,elevTable,nel,xt)
                newY(ncomp,j)= newY(1,linknb)

                xt = newY(ncomp,j)
                elevTable = xsec_tab(1,:,ncomp,j)
                areaTable = xsec_tab(2,:,ncomp,j)
                newArea(ncomp,j)=r_interpol(elevTable,areaTable,nel,xt)
                !dmy=area(1,linknb)+dac(1,linknb)

                dac(ncomp,j)=2*(newArea(ncomp,j)-oldArea(ncomp,j))-dap(ncomp,j)
                !* dqc(ncomp,j)
                dqc(ncomp,j)=dqc(1,linknb)*oldQ(ncomp,j)/oldQ(1,linknb)

                !* p.120,RM3
                qsum= 0.0
                linknb_ds= linknb
                do k=1, ndep(linknb_ds)
                    !* uslinks(k,j): k_th link ID that is immediately upstream of link j
                    linknb_us=uslinks(k,linknb_ds); nodenb=nx1(linknb_us)
                    qsum= qsum + qp(nodenb,linknb_us)
                end do
                qnp1_ds= oldQ(1,linknb_ds) +0.5*(dqp(1,linknb_ds)+dqc(1,linknb_ds))
                !* est. q(n+1, ncomp, link j_i), p120_RM
                qnp1_us= qnp1_ds*qp(ncomp,j)/qsum
                dqc(ncomp,j)= 2.0*(qnp1_us - oldQ(ncomp,j)) - dqp(ncomp,j)
            end if















        thes=thesinv
        call matrixc(j)


        do i=ncomp-1,1,-1
            !cour=dt(i)/dx(i)
            cour=dtini/dx(i,j)
            !rhs1=-cour*(f1(i+1)-f1(i)-d1(i+1)+d1(i))+lateralFlow(i)*dtini
            rhs1=-cour*(f1(i+1)-f1(i)-d1(i+1)+d1(i))+lateralFlow(i+1)*dtini
            rhs2=-cour*(f2(i+1)-f2(i)-d2(i+1)+d2(i))+dt(i)*grav*(ci2(i)+aso(i))
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



        if ( S_0 .gt. 0.) then
          if( y_crit_us .ge. y_norm_us ) then ! Steep slope
            frus2 = abs(newQ(1,j))/sqrt(grav*( dac(1,j)+areap(1,j) )**3.0/bo(1))
            print*, 'frus2', frus2, newQ(1,j), oldArea(1,j), areap(1,j), dac(1,j), bo(1)
            !pause 5000
            if (frus2 .ge. 1) then   ! Fr_us>=1: supercritical flow
              dqc(1,j)=dqp(1,j)
              dac(1,j)=dap(1,j)
            else
              dqc(1,j)=dqp(1,j)	!checked email
              dap(1,j)=dac(1,j)	!checked email
            endif
          else
            dqc(1,j)=dqp(1,j)	!checked email
            dap(1,j)=dac(1,j)	!checked email
          endif
        else
          dqc(1,j)=dqp(1,j)	!checked email
          dap(1,j)=dac(1,j)	!checked email
        endif


        !print*, newArea(35,1),area(35),areap(35,1),'65'

        ! Final update
        do i=1,ncomp
            da=(dap(i,j)+dac(i,j))/2.0
            dq=(dqp(i,j)+dqc(i,j))/2.0
            newArea(i,j)=da+oldArea(i,j)
            if(newArea(i,j) <= 0.0) newArea(i,j)=0.001

!           Now calculate y based on area calculated
!-------------------------------------
            elevTable = xsec_tab(1,:,i,j)
            areaTable = xsec_tab(2,:,i,j)

    !       interpolate the cross section attributes based on FINAL CALCULATED area
            xt=newArea(i,j)
            !print*, i, j, ncomp, newArea(35,1),'651'
            newY(i,j)=r_interpol(areaTable,elevTable,nel,xt)
            !print*, i, j, ncomp, newArea(35,1),'652'
!-------------------------------------

            newQ(i,j)=oldQ(i,j)+dq
            !froud(i)=abs(newQ(i))/sqrt(grav*newArea(i)**3.0/bo(i))

        end do
        !print*, newArea(1)
        !pause 5001

        !print*, newArea(35,1),'66'
        do i=1,ncomp
            currentQ = qp(i,j)
            q_sk_multi = 1.0
            do pp = 1, noQSKtable(j)
                if (  ( eachQSKtableNodeRange(1,pp,j) - i) * ( eachQSKtableNodeRange(2,pp,j) - i) .le. 0 ) then
                    tableLength = Q_sk_tableEntry(pp,j)
                    q_sk_multi = r_interpo_nn(Q_Sk_Table(1,1:tableLength,pp,j), &
                        Q_Sk_Table(2,1:tableLength,pp,j),tableLength,currentQ)
                end if
            end do
            co(i) = q_sk_multi*co(i)

            sfi = ( newQ(i,j) / co(i) ) ** 2.0
            celerity2(i)=5.0 / 3.0 * sfi ** 0.3 * abs(newQ(i,j)) ** 0.4 / bo(i) ** 0.4 / (1/(sk(i,j)*q_sk_multi)) ** 0.6
        end do
        celerity =  sum(celerity2) / ncomp


end subroutine mesh_dynamic_corrector
