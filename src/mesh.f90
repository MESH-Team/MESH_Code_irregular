!            FINITE DIFFERENCE METHOD
!
!  A program for one dimensional flow in open channel
!
program mesh

    use constants_module
    use arrays_module
    use arrays_section_module
    use var_module
    use matrix_module
    use sgate_module
    use xsec_attribute_module

    implicit none

    ! Local storage
    integer(kind=4) :: i, n, ntim, igate, pp, boundaryFileMaxEntry, ppp, qqq, saveFrequency !, areanew
    real(kind=4) :: cour, da, dq, dxini, yy, x, thetas, thesinv, tfin
    real(kind=4) :: t, skk, qq, qn, xt, r_interpol, t1, t2, maxCourant
    real(kind=4) :: arean, areac, hyrdn, hyrdc, perimn, perimc, qcrit, s0ds, timesDepth

    character(len=128) :: upstream_path , downstream_path
    character(len=128) :: bed_elevation_path, output_path, other_input
    character(len=128) :: channel_width_path
    character(len=128) :: path

    ! open file for input data
    character(len=128) :: dx_path

    call cpu_time( t1 )

    open(unit=1,file="../lower_Mississippi/input/input.txt",status='unknown')

    ! read data
    read(1,*) dtini
    read(1,*) dxini
    read(1,*) tfin
    read(1,*) ncomp
    read(1,*) ntim
    read(1,*) phi
    read(1,*) theta
    read(1,*) thetas
    read(1,*) thesinv
    read(1,*) alfa2
    read(1,*) alfa4
    read(1,*) f
    read(1,*) skk
    read(1,*) yy
    read(1,*) qq
    read(1,*) cfl
    read(1,*) ots
    read(1,*) yw
    read(1,*) bw
    read(1,*) w
    read(1,*) option
    read(1,*) yn
    read(1,*) qn
    read(1,*) igate
    read(1,*) xSection_path
    read(1,*) bed_elevation_path
    read(1,*) upstream_path
    read(1,*) downstream_path
    read(1,*) channel_width_path
    read(1,*) dx_path
    read(1,*) output_path
    read(1,*) option_dsbc
    read(1,*) maxTableLength
    read(1,*) nel
    read(1,*) timesDepth
    read(1,*) other_input
    read(1,*) boundaryFileMaxEntry
    read(1,*) saveFrequency
    close (1)

    ! Allocate arrays
    call setup_arrays(ntim, ncomp, maxTableLength, boundaryFileMaxEntry)
    call setup_arrays_section
    call setup_xsec_attribute_module(nel, ncomp)

    dt = dtini

    open(unit=90, file=trim(dx_path))
    do i=1,ncomp-1
        read(90, *) x, dx(i)
    end do
    close(90)

    ! reading Strickler's coefficient at each section
    open(unit=85,file=trim(other_input)//'Mannings_Stricklers_coeff.txt', status='unknown')
    do i=1,ncomp
        read(85, *) x, sk(i)
        call readXsection(i,(1.0/sk(i)),timesDepth)
        ! This subroutine creates attribute table for each cross sections and saves in the hdd
        ! setting initial condition
        !y(1,i) = yy ! + z(i)
        oldY(i) = yy ! + z(i)
    end do
    close(85)

    do i=1,ncomp
        call create_I2(i,ncomp)
    end do

    ityp = 1

    ! setting initial condition
    ! setting initial condition from previous work
    !open(unit=91,file=trim(output_path)//'initialCondition.txt', status='unknown')
    ! read(91, *)
    !do i=1,ncomp
    !    read(91, *) oldQ(i), oldY(i)
    !end do
    !close(91)

    !q(1, :) = qq
    oldQ = qq


    ! reading Q-Strickler's coefficient multiplier table
    open(unit=86,file=trim(other_input)//'Q_Mannings_table.txt', status='unknown')
    do i=1,maxTableLength
        read(86,*,end=300)Q_Sk_Table(1,i), Q_Sk_Table(2,i)
    enddo
300 close(86)
	Q_sk_tableEntry = i-1


    x = 0.0

    ! Read hydrograph input Upstream
    open(unit=87, file=upstream_path)
    do n=1,boundaryFileMaxEntry
        read(87,*,end=301) USBoundary(1, n), USBoundary(2, n)
    end do
301 close(87)
    ppp = n-1

    ! Read hydrograph input Downstream
    open(unit=88, file=downstream_path)
    do n=1,boundaryFileMaxEntry
      read(88,*,end=302)  DSBoundary(1, n), DSBoundary(2, n)
    end do
302 close(88)
    qqq = n-1

    t = 0.0

    ! applying boundary
    ! interpolation of boundaries at the initial time step

    oldQ(1)    =r_interpol(USBoundary(1, :),USBoundary(2, :),ppp,t)
    oldY(ncomp)=r_interpol(DSBoundary(1, :),DSBoundary(2, :),qqq,t)

    ! DS Boundary treatment: from water level to area time series
    ncompElevTable = xsec_tab(1,:,ncomp)
    ncompAreaTable = xsec_tab(2,:,ncomp)

	!open(unit=81,file=trim(output_path)//'DS_area.txt', status='unknown')
    xt=oldY(ncomp)
    oldArea(ncomp)=r_interpol(ncompElevTable,ncompAreaTable,nel,xt)
    !write(81, *) t, oldArea(ncomp)

    ! Open files for output
    path = trim(output_path) // 'output_wl.txt'
    open(unit=8, file=trim(path), status='unknown')
    path = trim(output_path) // 'q.txt'
    open(unit=9, file=trim(path), status='unknown')

    path = trim(output_path) // 'area.txt'
    open(unit=51, file=trim(path), status='unknown')

    ! Output initial condition

    write(8, 10)  t, (oldY(i), i=1,ncomp)
    write(9, 10)  t, (oldQ(i), i=1, ncomp)
    write(51, 10) t, (oldArea(i), i=1, ncomp)

    !
    ! Loop on time
    !
    do n=1, ntim-1

        ! interpolation of boundaries at the desired time step
        newQ(1)     =r_interpol(USBoundary(1, :),USBoundary(2, :),ppp,t+dtini)
        newY(ncomp) =r_interpol(DSBoundary(1, :),DSBoundary(2, :),qqq,t+dtini)

        xt=newY(ncomp)
		newArea(ncomp)=r_interpol(ncompElevTable,ncompAreaTable,nel,xt)

        ! Set upstream discharge
        dqp(1) = newQ(1) - oldQ(1)

        call section()
        ! Nazmul: The subroutine calls the attribute tables and interpolate according to the available water level
        thes=thetas

        call matrixp()

        do i=2,ncomp
            cour=dt(i)/dx(i-1)
            rhs1=-cour*(f1(i)-f1(i-1)-d1(i)+d1(i-1))
            rhs2=-cour*(f2(i)-f2(i-1)-d2(i)+d2(i-1))+dt(i)*grav*(ci2(i)+aso(i))
            c11=g11inv(i)*b11(i-1)+g12inv(i)*b21(i-1)
            c12=g11inv(i)*b12(i-1)+g12inv(i)*b22(i-1)
            c21=g21inv(i)*b11(i-1)+g22inv(i)*b21(i-1)
            c22=g21inv(i)*b12(i-1)+g22inv(i)*b22(i-1)
            dap(i)=g11inv(i)*rhs1+g12inv(i)*rhs2-c11*dap(i-1)-c12*dqp(i-1)
            dqp(i)=g21inv(i)*rhs1+g22inv(i)*rhs2-c21*dap(i-1)-c22*dqp(i-1)
        end do

        ! Boundary conditions at downstream (right boundary)
        if (option_dsbc.eq.1) then
            dac(ncomp)=dap(ncomp)
            yn=(area(ncomp)+dap(ncomp))/bo(ncomp)
            arean=yn*bo(ncomp)
            perimn=2.0*yn+bo(ncomp)
            hyrdn=arean/perimn
            s0ds=-((z(ncomp)-z(ncomp-1))/dx(ncomp))
            qn=skk*arean*hyrdn**(2.0/3.0)*sqrt(s0ds)
            dqp(ncomp)=qn-oldQ(ncomp)
            dqc(ncomp)=dqp(ncomp)

        elseif(option_dsbc.eq.2)then
            dac(ncomp)=dap(ncomp)
            yn=(area(ncomp)+dap(ncomp))/bo(ncomp)
            areac=yn*bo(ncomp)
            perimc=2.0*yn+bo(ncomp)
            hyrdc=areac/perimc
            s0ds=-((z(ncomp)-z(ncomp-1))/dx(ncomp))
            qn=skk*areac*hyrdc**(2.0/3.0)*sqrt(s0ds)
            qcrit=1.05*(((yn**3.0)*(bo(ncomp)**2.0)*grav)**(1.0/2.0))
            write(*,*)qcrit
            dqp(ncomp)=qcrit-oldQ(ncomp)
            dqc(ncomp)=dqp(ncomp)

        else
            !dac(ncomp)=0.0
            !dap(ncomp)=0.0
            !dqc(ncomp)=dqp(ncomp)

!            dap(ncomp)=0.0	!checked email !for critical comment out
! change for unsteady flow
			dap(ncomp) = newArea(ncomp) - oldArea(ncomp)
            dac(ncomp)=dap(ncomp)	!checked email
            dqc(ncomp)=dqp(ncomp)	!checked email

        endif

        ! Update via predictor
        areap = area + dap
        qp = oldQ + dqp

        call secpred()
        thes=thesinv
        call matrixc()

        do i=ncomp-1,1,-1
            cour=dt(i)/dx(i)
            rhs1=-cour*(f1(i+1)-f1(i)-d1(i+1)+d1(i))
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

        dqc(1)=dqp(1)	!checked email
        !dac(1)=dap(1) !for critical, uncomment
        dap(1)=dac(1)	!checked email

        ! Final update
        do i=1,ncomp
            da=(dap(i)+dac(i))/2.0
            dq=(dqp(i)+dqc(i))/2.0
            newArea(i)=da+area(i)
            if(newArea(i) <= 0.0) newArea(i)=0.001

!           Now calculate y based on area calculated
!-------------------------------------
            do pp=1,nel
                elevTable(pp) = xsec_tab(1,pp,i)
                areaTable(pp) = xsec_tab(2,pp,i)
            enddo

    !       interpolate the cross section attributes based on FINAL CALCULATED area
            xt=newArea(i)
            newY(i)=r_interpol(areaTable,elevTable,nel,xt)
!-------------------------------------

            newQ(i)=oldQ(i)+dq
            froud(i)=abs(newQ(i))/sqrt(grav*newArea(i)**3.0/bo(i))

        end do

        do i=1,ncomp-1
            courant(i)=(newQ(i)+newQ(i+1))/(newArea(i)+newArea(i+1))*dtini/dx(i)
        enddo

        if (maxCourant .lt. maxval (courant)) then
            maxCourant = maxval (courant)
        endif

        t = t + dtini
        print "('- cycle',i9,'  completed')", n

        if (mod(n,saveFrequency) .eq. 0) then
        write(8, 10) t, (newY(i), i=1,ncomp)
        write(9, 10) t, (newQ(i), i=1,ncomp)
        write(51, 10) t, (newArea(i), i=1, ncomp)
        end if

        !write(81, *) t, newArea(ncomp)

        ! update of Y, Q and Area vectors
        oldY   = newY
        oldQ   = newQ
        oldArea= newArea

    end do
    ! End of time loop

    close(8)
    close(9)
    close(51)
    !close(81)

    print*, 'dx', (dx(i), i=1, ncomp-1)
    print*, 'Froude', (froud(i), i=1, ncomp)
    print*, 'Bed', (z(i), i=1, ncomp)
    print*, 'newArea', (newArea(i), i=1, ncomp)
    print*, 'I2_corr', (ci2(i), i=1, ncomp)
    print*, 'Courant no', (courant(i), i=1, ncomp-1)
    print*, 'Maximum Courant no', maxCourant

    !
10  format(f12.2 , <ncomp>f12.2)
!10  format(f12.2 , 1000f12.2)

    call cpu_time( t2 )
    print '("Time = ",f10.3," seconds.")',t2 - t1
    pause 202
end program mesh