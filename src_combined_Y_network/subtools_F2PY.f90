module subtools

    use constants_module
    use arrays_module
    use var_module
    use arrays_section_module
    use xsec_attribute_module

contains


    subroutine calc_q_sk_multi(i,j,currentQ,q_sk_multi)
        implicit none

        integer,intent(in) :: i, j
        real,intent(in) :: currentQ
        real,intent(out) :: q_sk_multi
        integer :: pp, tableLength
        real :: r_interpo_nn

        q_sk_multi = 1.0
        do pp = 1, noQSKtable(j)
            if (  ( eachQSKtableNodeRange(1,pp,j) - i) * ( eachQSKtableNodeRange(2,pp,j) - i) .le. 0 ) then
                tableLength = Q_sk_tableEntry(pp,j)
                q_sk_multi = r_interpo_nn(Q_Sk_Table(1,1:tableLength,pp,j), &
                Q_Sk_Table(2,1:tableLength,pp,j),tableLength,currentQ)
            end if
        end do

    end subroutine


    subroutine r_interpol(x,y,kk,xrt,yt)
        implicit none

    integer, intent(in) :: kk
    real, intent(in) :: xrt, x(kk), y(kk)
    real, intent(out) :: yt
    integer :: k

    if (xrt.le.maxval(x) .and. xrt.ge.minval(x)) then
        do k=1,kk-1
            if((x(k)-xrt)*(x(k+1)-xrt).le.0)then

                yt=(xrt-x(k))/(x(k+1)-x(k))*(y(k+1)-y(k))+y(k)

                EXIT
            endif
        end do
    else
        print*, xrt, ' is not within the limit'
        yt = -9999.0
        print*, 'maxval(x)= ', maxval(x), 'and minval(x)=', minval(x),'so',  xrt, ' is not within the limit'
        print*, 'kk', kk
        print*, 'x', (x(k), k=1, kk)
        print*, 'y', (y(k), k=1, kk)
        !stop
        !if (xrt.le. minval(x)) yt=minval(y)
        !if (xrt.ge. maxval(x)) yt=maxval(y)
    end if
    !r_interpol = yt
    !return
end subroutine r_interpol


!+++-------------------------------------------------------------------
!+ Compute Q from a given value of Y using Manning's equation, which
!+ becomes uniform flow eqns when slp is So and nonuniform flow eqns
!+ when slp is Sf.
!+++-------------------------------------------------------------------
    subroutine ManningEq_QforY(y_mn, bo_mn, slp, n_mn, q_mn)

        implicit none

        real,intent(in) :: y_mn, bo_mn, slp, n_mn
        real,intent(out) :: q_mn
        real :: ar, peri, hydr
        integer :: chshp

        !+++----------------------------------------------------+
        !+ Rectangular channels
        !+++----------------------------------------------------+
        if (chshp.eq.1) then
            !* area
            ar=y_mn*bo_mn
            !* perimeter
            peri=2.0*y_mn+bo_mn
            !* hydraulic radius
            hydr=ar/peri
            !* Note that n_mn=1/Mannin's N.
            !q_mn=n_mn*ar*(hydr**(2.0/3.0))*sqrt(slp)
            q_mn=n_mn*ar*(hydr**(2.0/3.0))*(slp**0.5)
        end if


    end subroutine ManningEq_QforY

!+++-------------------------------------------------------------------------------
!+ Compute Y and thus area from a given value of Q using Manning's equation, which
!+ becomes uniform flow eqns when slp is So and nonuniform flow eqns
!+ when slp is Sf. The numerical method used in this subroutine is described in
!+ p.88, Chaudhry book.
!+++-------------------------------------------------------------------------------
    !subroutine ManningEq_YforQ()!(y_mn, bo_mn, slp, n_mn, q_mn)

    !end subroutine ManningEq_YforQ

!+++-------------------------------------------------------------------------------
!+ computation of normal depth in regular/trapezoidal x-section using
!+ Newton-Raphson method. Refer to Appendix C-2, Chaudhary and p71,RM1_MESH
!+++-------------------------------------------------------------------------------
    !subroutine normal_crit_y(j,ynm0, q_sk_multi, So, dsc, y_norm, y_crit, area_n, area_c)
    subroutine normal_crit_y(i, j, q_sk_multi, So, dsc, y_norm, y_crit, area_n, area_c)


        implicit none

        integer, intent(in) :: i, j
        real, intent(in) :: q_sk_multi, So, dsc
        real, intent(out) :: y_norm, y_crit, area_n, area_c
        real :: area_0, width_0, errorY, hydR_0!, fro
        integer :: trapnm_app, recnm_app, iter


            elevTable = xsec_tab(1,:,i,j)
            areaTable = xsec_tab(2,:,i,j)
            rediTable = xsec_tab(4,:,i,j)
            topwTable = xsec_tab(6,:,i,j)
            !print*, 'initial ara 0', oldY(i,j)
            call r_interpol(elevTable,areaTable,nel,oldY(i,j),area_0) ! initial estimate
            if (area_0 .eq. -9999) then
                print*, 'At j = ',j,', i = ',i, 'interpolation of area_0 in calculating normal area was not possible'
                stop
            end if
            call r_interpol(elevTable,topwTable,nel,oldY(i,j),width_0) ! initial estimate

            !print*, 'initial ara 0', area_0

            area_c=area_0
            errorY = 100.
            !pause
            do while (errorY .gt. 0.00001)
                !print*, '11', area_0
                call r_interpol(areaTable,rediTable,nel,area_0,hydR_0)
                if (hydR_0 .eq. -9999) then
                    print*, 'At j = ',j,', i = ',i, 'interpolation of hydR_0 in calculating normal area was not possible'
                    stop
                end if
                !print*, '12', hydR_0
                area_n = dsc/sk(i,j)/q_sk_multi/ hydR_0 ** (2./3.) / sqrt(So)
                !print*, '13', dsc, area_n

                errorY = abs(area_n - area_0) / area_n
                area_0 = area_n
                !print*, 'area_0', area_n,'hydR_0', hydR_0
                area_c = (dsc * dsc * width_0 / grav) ** (1./3.)
                call r_interpol(areaTable,topwTable,nel,area_c, width_0)
                !fro=abs(dsc)/sqrt(grav*area_c**3.0/width_0)
            enddo

            call r_interpol(areaTable,elevTable,nel,area_0,y_norm)
            if (y_norm .eq. -9999) then
                print*, 'At j = ',j,', i = ',i, 'interpolation of y_norm in calculating normal area was not possible'
                stop
            end if
            call r_interpol(areaTable,elevTable,nel,area_c,y_crit)

    end subroutine normal_crit_y

    subroutine normal_y(i, j, q_sk_multi, So, dsc, y_norm, y_crit, area_n, area_c)


        implicit none

        integer, intent(in) :: i, j
        real, intent(in) :: q_sk_multi, So, dsc
        real, intent(out) :: y_norm, y_crit, area_n, area_c
        real :: area_0, width_0, errorY, hydR_0!, fro
        integer :: trapnm_app, recnm_app, iter


            elevTable = xsec_tab(1,:,i,j)
            areaTable = xsec_tab(2,:,i,j)
            rediTable = xsec_tab(4,:,i,j)
            topwTable = xsec_tab(6,:,i,j)
            !print*, 'initial ara 0', oldY(i,j)
            call r_interpol(elevTable,areaTable,nel,oldY(i,j),area_0) ! initial estimate
            !call r_interpol(elevTable,topwTable,nel,oldY(i,j),width_0) ! initial estimate
            !print*, 'here2'

            !print*, 'initial ara 0', area_0

            !area_c=area_0
            errorY = 100.
            !pause
            do while (errorY .gt. 0.00001)
                !print*, '11', area_0
                call r_interpol(areaTable,rediTable,nel,area_0,hydR_0)
                !print*, '12', hydR_0
                area_n = dsc/sk(i,j)/q_sk_multi/ hydR_0 ** (2./3.) / sqrt(So)
                !print*, '13', dsc, area_n

                errorY = abs(area_n - area_0) / area_n
                area_0 = area_n
                !print*, 'area_0', area_n,'hydR_0', hydR_0
                !area_c = (dsc * dsc * width_0 / grav) ** (1./3.)
                !call r_interpol(areaTable,topwTable,nel,area_c, width_0)

                !fro=abs(dsc)/sqrt(grav*area_c**3.0/width_0)
            enddo
            call r_interpol(areaTable,elevTable,nel,area_0,y_norm)
            !call r_interpol(areaTable,elevTable,nel,area_c,y_crit)
            y_crit = -9999.
            area_c = -9999.

    end subroutine normal_y

!+++----------------------------------------------------------------------------
!+ computation of critical depth in regular/trapezoidal x-section using
!+ Newton-Raphson method. Refer to Appendix B-2, Chaudhary and p69-70,RM1_MESH
!+++----------------------------------------------------------------------------
    subroutine criticaldep_NR(sslp, bwd, dsc, ycr)

        implicit none

        real, intent(in) :: sslp, bwd, dsc
        real, intent(out) :: ycr
        real :: Axs, Bxs, c1, dbdy, fy, dfdy, ycrnew
        real :: errNR, tolNR
        real :: dmy, epsc, tc0, etac
        integer :: trapcr_app, chshp

        tolNR=0.01  !*tolerance
        errNR=10.0*tolNR

        trapcr_app =2   !*approach 1 or 2 for trapezoidal channel


        if (chshp==1) then
        !* rectangular channel
            ycr=(dsc**(2.0/3.0))/((grav*(bwd**2.0))**(1.0/3.0))

        elseif (chshp==2) then
        !* trapezoidal channel
            if (trapcr_app == 1) then
            !* approach 1 (Newton-Rapson)
                ycr=(dsc**(2.0/3.0))/((grav*(bwd**2.0))**(1.0/3.0)) !*initial estimate
                !c1=dsc/sqrt(grav)
                c1=dsc/(grav**0.5)
                dbdy=2.0*sslp

                do while (errNR>tolNR)
                    Axs=(bwd+sslp*ycr)*ycr
                    Bxs=bwd+2.0*sslp*ycr
                    !fy=sqrt((Axs**3.0)/Bxs)-c1
                    fy=((Axs**3.0)/Bxs)**0.5 - c1

                    !dfdy=1.5*sqrt(Axs*Bxs)-0.5*((Axs/Bxs)**1.5)*dbdy
                    dfdy=1.5*((Axs*Bxs)**0.5)-0.5*((Axs/Bxs)**1.5)*dbdy
                    ycrnew=ycr-fy/dfdy
                    errNR = abs((ycrnew-ycr)/ycrnew)
                    ycr=ycrnew
                end do
            !* approach 2 (based on "Explicit solutions for critical and normal depths
            !* in trapezoidal and parabolic open channels" by Ali R.V.
            elseif (trapcr_app == 2) then
                epsc=4.0*(((sslp**3.0)*(dsc**2.0)/(grav*(bwd**5.0)))**(1.0/3.0))

                dmy=(1.0+0.666*(epsc**1.041))**0.374
                tc0=(1.0+1.161*epsc*dmy)**0.144

                dmy=(5.0*(tc0**6.0) + 1.0)/(6.0*(tc0**5.0) - epsc)
                !dmy=(5.0*tc0**0.864 + 1)/(6.0*tc0**0.72 - epsc)
                etac=-0.5+0.5*(dmy**3.0)
                ycr=bwd*etac/sslp
            end if
        end if

    end subroutine criticaldep_NR


!+++----------------------------------------------------------------------------
!+ computation of downstream conjugate depth for regular/trapezoidal x-section
!+ Refer to p.40,Chaudhary and p72,RM1_MESH
!+++----------------------------------------------------------------------------
    subroutine conjugatedep(ycj1 ,sslp, bwd, dsc, ycj2)
        implicit none
        real, intent(in) :: ycj1, sslp, bwd, dsc
        real, intent(out) :: ycj2
        real :: Axs, Fr1
        real uwQ, xcj, eta, kcj, delta, tcj, ycj0, lamX
        integer iter, mxiter, trapcj_app
        integer :: chshp

        trapcj_app=2    !*approach 1 or 2 for trapezoidal channel

        if (chshp==1) then
        !* rectangular channel
            Axs=ycj1*bwd
            !Fr1=dsc/sqrt(grav*(Axs**3.0)/bwd)
            Fr1=dsc/((grav*(Axs**3.0)/bwd)**0.5)
            !ycj2=0.5*ycj1*(-1.0+sqrt(1.0+8.0*(Fr1**2.0)))
            ycj2=0.5*ycj1*(-1.0 + (1.0+8.0*(Fr1**2.0))**0.5)
        elseif (chshp==2) then
        !* trapezoidal channel
            !* Two approaches based on "Computing conjugate depths in trapezoidal channels" by Liu J. et al (2012).
            !* For more details, refer to p73~74, RM1_MESH
            uwQ=dsc/bwd                         !*unit width discharge
            xcj=ycj1/(uwQ**(2.0/3.0))
            eta=sslp*(uwQ**(2.0/3.0))/bwd         !*eta
            kcj=4.0*eta*((1.0/grav)**(1.0/3.0))

            delta=(((1.0+kcj*((1.0+kcj)**0.2))**0.5)-1.0)/2.0
            tcj=delta/eta

            ycj0=tcj + (tcj-xcj)*((tcj/xcj)**0.5)
            lamX=6.0/(grav*xcj*(1+eta*xcj))+(xcj**2.0)*(3.0+2.0*eta*xcj)

            if (trapcj_app==1) then
            !* Approach 1 (iterative formula)
                mxiter=3
                do iter=1,mxiter
                    ycj0=((lamX-6.0/(grav*ycj0*(1.0+eta*ycj0)))/(3.0+2.0*eta*ycj0))**0.5
                end do
            elseif (trapcj_app==2) then
            !* Approach 2 (estimation formula)
                ycj0=((lamX-6.0/(grav*ycj0*(1.0+eta*ycj0)))/(3.0+2.0*eta*ycj0))**0.5
            end if

            ycj2=ycj0*(uwQ**(2.0/3.0))
        end if

    end subroutine conjugatedep

!+++-----------------------------------------------------
!+ Computation of depth with given area, p77, RM1_MESH
!+++-----------------------------------------------------
    subroutine depthcalc(sslp, bwd, Axs, dpth)

        implicit none

        real, intent(in) :: sslp, bwd, Axs
        real, intent(out) :: dpth
        integer :: chshp

        if (chshp==1) then
        !* rectangular channel
            dpth=Axs/bwd
        elseif (chshp==2) then
        !* trapezoidal channel
            !dpth=(-bwd + sqrt((bwd**2.0) + 4.0*sslp*Axs))/(2.0*sslp)
            dpth=(-bwd + ((bwd**2.0) + 4.0*sslp*Axs)**0.5)/(2.0*sslp)
        end if

    end subroutine depthcalc

!+++-----------------------------------------------------
!+ Computation of area of various channel x-sections
!+++-----------------------------------------------------
    subroutine areacalc(yxs, sslp, bwd, Axs)
          implicit none

        real, intent(in) :: yxs, sslp, bwd
        real, intent(out) :: Axs
        integer :: chshp

        if (chshp==1) then
        !* rectangular channel
            Axs=bwd*yxs
        elseif (chshp==2) then
        !* trapezoidal channel
            Axs=(bwd+sslp*yxs)*yxs
        end if

    end subroutine areacalc

!+++-----------------------------------------------------
!+ Computation of I1 (=centroid*A), p80,RM1_MESH
!+++-----------------------------------------------------
    subroutine I1calc(yxs, sslp, bwd, Axs, I1)

        implicit none

        real, intent(in) :: yxs, sslp, bwd, Axs
        real, intent(out) :: I1
        integer :: chshp

        if (chshp==1) then
        !* rectangular channel
            I1=0.5*bwd*(yxs**2.0)
        elseif (chshp==2) then
        !* trapezoidal channel
            I1=yxs*(3.0*bwd + 4.0*sslp*yxs)*Axs/(6.0*(bwd+sslp*yxs))
        end if

    end subroutine I1calc


!+++-----------------------------------------------------
!+ Computation of hydraulic radius R (=A/P)
!+++-----------------------------------------------------
    subroutine hydRcalc(yxs, sslp, bwd, Axs, hydR)

        implicit none

        real, intent(in) :: yxs, sslp, bwd, Axs
        real, intent(out) :: hydR
        integer :: chshp, gdI2dA

        if (chshp==1) then
        !* rectangular channel
            hydR=Axs/(bwd+2.0*yxs)
        elseif (chshp==2) then
        !* trapezoidal channel
            !hydR=Axs/(bwd+2.0*yxs*sqrt(1.0+(sslp**2.0)))
            hydR=Axs/(bwd+2.0*yxs*((1.0+(sslp**2.0))**0.5))
        end if

    end subroutine hydRcalc

!+++-----------------------------------------------------
!+ Computation of c in matrix A,p80~80-2,RM1_MESH
!+++-----------------------------------------------------
    subroutine c_mtrxAcalc(yxs, sslp, bwd, Axs, cA)

        implicit none

        real, intent(in) :: yxs, sslp, bwd, Axs
        real, intent(out) :: cA
        real :: b24sa, h, sxbsh, dhdA, the1, the2, dycdA, yc, dycAdA
        integer :: chshp, gdI2dA


        if (chshp==1) then
        !* rectangular channel
            !cA=sqrt(grav*Axs/bwd)
            cA=(grav*Axs/bwd)**0.5
        elseif (chshp==2) then
        !* trapezoidal channel
            b24sa=(bwd**2.0) + 4.0*sslp*Axs
            h=yxs !(-bwd + sqrt(b24sa))/(2.0*sslp)
            sxbsh=6.0*(bwd + sslp*h)
            dhdA= (b24sa)**(-0.5)
            the1=(3.0*bwd + 8.0*sslp*h)*dhdA
            the2=6.0*sslp*dhdA
            dycdA=(the1-h*the2*(3.0*bwd+4.0*sslp*h)/sxbsh)/sxbsh

            yc=h*(3.0*bwd+4.0*sslp*h)/(6.0*(bwd+sslp*h))
            dycAdA=(yc + Axs*dycdA)

            !cA=sqrt(grav*dycAdA)
            cA=(grav*dycAdA)**0.5
        end if

    end subroutine c_mtrxAcalc


!+++-----------------------------------------------------
!+ Computation of dk/dA,p82-83,RM1_MESH
!+++-----------------------------------------------------
    subroutine dKdAcalc(manN, sslp, bwd, Axs, dkda)

        implicit none

        real, intent(in) :: manN, sslp, bwd, Axs
        real, intent(out) :: dkda
        real :: eta, pi, phi
        integer :: chshp, gdI2dA

        if (chshp==1) then
        !* rectangular channel
            eta=(bwd+2.0*Axs/bwd)
            dkda=((5.0/3.0)*(Axs**(2.0/3.0))*eta - (4.0/3.0)*(Axs**(5.0/3.0))/bwd)/(manN*(eta**(5.0/3.0)))
            !dkda=((5.0/3.0)*(Axs**(2.0/3.0))*eta - 2.0*(Axs**(5.0/3.0))/bwd)/(manN*eta**2.0)
        elseif (chshp==2) then
        !* trapezoidal channel
            !pi=bwd + (-bwd + sqrt((bwd**2.0) + 4.0*sslp*Axs))*sqrt(1.0+(sslp**2.0))/sslp
            pi=bwd + (-bwd + ((bwd**2.0) + 4.0*sslp*Axs)**0.5)*((1.0+(sslp**2.0))**0.5)/sslp
            !phi=2.0*sqrt(1.0+(sslp**2.0))*(((bwd**2.0) + 4.0*sslp*Axs)**(-0.5))
            phi=2.0*((1.0+(sslp**2.0))**0.5)*(((bwd**2.0) + 4.0*sslp*Axs)**(-0.5))
            dkda=((5.0/3.0)*(Axs**(2.0/3.0))*pi - (2.0/3.0)*(Axs**(5.0/3.0))*phi)/(manN*(pi**(5.0/3.0)))
        end if

    end subroutine dKdAcalc


!+++-----------------------------------------------------
!+ Computation of g*dI2/dA,p81,RM1_MESH
!+++-----------------------------------------------------
    subroutine gdI2dAcalc(sslp, bwd, Axs, dsgmdx, gdI2dA)

        implicit none

        real, intent(in) :: sslp, bwd, Axs, dsgmdx
        real, intent(out) :: gdI2dA
        integer :: chshp

        if (chshp==1) then
        !* rectangular channel
            gdI2dA=grav*(Axs/(bwd**2.0))*dsgmdx
        elseif (chshp==2) then
        !* trapezoidal channel
            gdI2dA=grav*(1.0-bwd*(((bwd**2.0) + 4.0*sslp*Axs)**(-0.5)))/(2.0*sslp)*dsgmdx
        end if

    end subroutine gdI2dAcalc


 end module
