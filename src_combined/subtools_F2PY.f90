NOT USED NOW
module subtools

    use constants_module
    use arrays_module
    use var_module
    use arrays_section_module
    use xsec_attribute_module

contains


 !   subroutine calc_q_sk_multi(i,currentQ,q_sk_multi)
 !       implicit none
 !
 !       integer,intent(in) :: i
 !       real,intent(in) :: currentQ
 !       real,intent(out) :: q_sk_multi
 !       integer :: pp, tableLength
 !       real :: r_interpo_nn

 !       q_sk_multi = 1.0
 !       do pp = 1, noQSKtable
 !           if (  ( eachQSKtableNodeRange(1,pp) - i) * ( eachQSKtableNodeRange(2,pp) - i) .le. 0 ) then
 !               tableLength = Q_sk_tableEntry(pp)
 !               q_sk_multi = r_interpo_nn(Q_Sk_Table(1,1:tableLength,pp), &
 !               Q_Sk_Table(2,1:tableLength,pp),tableLength,currentQ)
 !           end if
 !       end do

 !   end subroutine


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
    subroutine normal_crit_y(i, q_sk_multi, So, dsc, y_norm, y_crit, area_n, area_c)


        implicit none

        integer, intent(in) :: i
        real, intent(in) :: q_sk_multi, So, dsc
        real, intent(out) :: y_norm, y_crit, area_n, area_c
        real :: area_0, width_0, errorY, hydR_0, r_interpol!, fro
        integer :: trapnm_app, recnm_app, iter


        !c1=manN*dsc/(c0*sqrt(So))

        !y_norm=1.0    !*initial estimate


            elevTable = xsec_tab(1,:,i)
            areaTable = xsec_tab(2,:,i)
            rediTable = xsec_tab(4,:,i)
            topwTable = xsec_tab(6,:,i)
            area_0 = r_interpol(elevTable,areaTable,nel,oldY(i))
            width_0= r_interpol(elevTable,topwTable,nel,oldY(i))

            area_c=area_0
            errorY = 100.
            !print*, 'inside function', i, q_sk_multi, So, dsc, area_0, width_0
            !pause
            do while (errorY .gt. 0.00001)

                hydR_0 = r_interpol(areaTable,rediTable,nel,area_0)
                area_n = dsc/sk(i)/q_sk_multi/ hydR_0 ** (2./3.) / sqrt(So)

                errorY = abs(area_n - area_0) / area_n
                area_0 = area_n
                area_c = (dsc * dsc * width_0 / grav) ** (1./3.)
                width_0  = r_interpol(areaTable,topwTable,nel,area_c)
                !fro=abs(dsc)/sqrt(grav*area_c**3.0/width_0)
                !print*, 'check point -1', area_0,  area_c, fro

            enddo
                !print*, 'nazmul j', j, oldY(1,j)

            y_norm = r_interpol(areaTable,elevTable,nel,area_0)
            y_crit = r_interpol(areaTable,elevTable,nel,area_c)


            !print*, 'check point 0', So, S_0,dsc,width_0, y_norm-z(i), y_crit-z(i), area_n, area_norm, area_c, fro
            !print*, 'check point -1', dsc, y_norm-z(1),  y_crit-z(1)
            !pause 500
    end subroutine normal_crit_y

 end module
