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
    integer(kind=4) :: i, j, n, ntim, igate, pp, boundaryFileMaxEntry, ppp, qqq, noLatFlow, noQSKtable
    real(kind=4) :: cour, da, dq, dxini, yy, x, thetas, thesinv
    real(kind=4) :: skk, qq, qn, xt, r_interpol, t1, t2, maxCourant, r_interpo_nn
    real(kind=4) :: arean, areac, hyrdn, hyrdc, perimn, perimc, qcrit, s0ds, timesDepth, latFlowValue
	real(kind=4) :: t, r_interpol_time, t0, saveInterval, tfin

    real(kind=4) :: q_sk_multi, sfi, soi

    character(len=128) :: upstream_path , downstream_path
    character(len=128) :: manning_strickler_path, output_path, other_input
    character(len=128) :: QSKtablePath, lateralFlow_path
    character(len=128) :: path

    !!! TEST !!!
    real(kind=4) :: a1, a2, a3, a4, b1, b2, b3, b4, dd1, dd2, dd3, dd4, h1, h2, h3, h4, alpha
    real(kind=4) :: qy, qxy, qxxy, qxxxy, ppi, qqi, rri, ssi, sxi, cour2, mannings, Sb, width

    real(kind=4) :: pere(500)
    !real(kind=4) ::



    ! open file for input data
    character(len=128) :: dx_path

    call cpu_time( t1 )

    !open(unit=1,file="../lower_Mississippi/input/input_BR_2_BC_2009.txt",status='unknown')
    !open(unit=1,file="../lower_Mississippi/input_BC2SWP/input_BC_2_HOP_2009.txt",status='unknown')
    !open(unit=1,file="../lower_Mississippi/input/input_test.txt",status='unknown')
    !open(unit=1,file="../Crank_Nicolson_project_BR2BC/input/input_BR_2_BC_2009.txt",status='unknown')
    !open(unit=1,file="../Crank_Nicolson_project/input/input_crank_nicolson.txt",status='unknown')
    !open(unit=1,file="../Crank_Nicolson_project/input/input_crank_nicolson_test_lateralFlow.txt",status='unknown')

    !open(unit=1,file="../Mississippi_River_11_years_20200511/input/input_Mississippi_BR2SWP_diffusive.txt",status='unknown')
    !open(unit=1,file="../Mississippi_River_11_years_20200511/input/input_Mississippi_BR2SWP_muskingumCunge.txt",status='unknown')
    open(unit=1,file="../Vermelion_River/input/input_Vermelion_diffusive_20200526.txt",status='unknown')
    !open(unit=1,file="../Vermelion_River/input/input_Vermelion_MuskingumCunge.txt",status='unknown')


    print*, 'Reading input file'

    ! read data
    read(1,*) dtini    !! in seconds
    read(1,*) dxini
    read(1,*) t0    !! in hours
    read(1,*) tfin    !! in hours
    read(1,*) ncomp
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
    read(1,*) manning_strickler_path
    read(1,*) upstream_path
    read(1,*) downstream_path
    read(1,*) QSKtablePath
    read(1,*) dx_path
    read(1,*) lateralFlow_path
    read(1,*) output_path
    read(1,*) option_dsbc
    read(1,*) maxTableLength
    read(1,*) nel
    read(1,*) timesDepth
    read(1,*) other_input
    read(1,*) boundaryFileMaxEntry
    read(1,*) saveInterval
    read(1,*) noLatFlow

    allocate(latFlowLocations(noLatFlow))   ! all the first nodes where a lateral flow starts
    allocate(latFlowType(noLatFlow))        ! Lateral flow type: Type = 1 = time series; Type 2 = flow as a function of upstream flow
    allocate(latFlowXsecs(noLatFlow))       ! no of x-secs at the downstream that the lateral flow is applied

    read(1,*) (latFlowLocations(i), i=1, noLatFlow)
    do i=1,noLatFlow
        if ((latFlowLocations(i)-1)*(latFlowLocations(i)-ncomp) .eq. 0) then
            print*, 'ERROR: Lateral flow cannot be applied at the boundaries'
            stop
        end if
    end do
    read(1,*) (latFlowType(i), i=1, noLatFlow)
    do i=1,noLatFlow
        if (latFlowType(i) .eq. 1) then
            print*, 'Lateral flow at node = ', latFlowLocations(i), ', is a time series'
        elseif (latFlowType(i) .eq. 2) then
            print*, 'Lateral flow at node = ', latFlowLocations(i), ', is a function of upstream flow'
        else
            print*, 'Wrong lateral flow type is provided. Type ', latFlowType(i), 'is not a valid type'
            stop
        end if
    end do

    read(1,*) (latFlowXsecs(i), i=1, noLatFlow)

    read(1,*) noQSKtable

    !allocate(QSKTableData(2,boundaryFileMaxEntry,noQSKtable))   ! all the first nodes where a lateral flow starts
    allocate(eachQSKtableNodeRange(2,noQSKtable))        ! Lateral flow type: Type = 1 = time series; Type 2 = flow as a function of upstream flow

    read(1,*) (eachQSKtableNodeRange(1,i), i=1, noQSKtable)
    read(1,*) (eachQSKtableNodeRange(2,i), i=1, noQSKtable)

    !!! Need to test so that one section does not corresponds to more than one table
    do i = 2, noQSKtable
        if ( eachQSKtableNodeRange(2,i-1) .ge. eachQSKtableNodeRange(1,i) ) then
            print*, 'Wrong range of nodes applied for Q-Sk table.'
            print*, 'Lower limit of Table ', i-1,'must be smaller than the upper limit of Table ', i
            stop
        end if
    end do


    close (1)
    ! reading input is completed


    ! Allocate arrays
    call setup_arrays(ntim, ncomp, maxTableLength, boundaryFileMaxEntry, noLatFlow, noQSKtable)
    call setup_arrays_section
    call setup_xsec_attribute_module(nel, ncomp)


    ! setting up the time steps
    dt = dtini


    open(unit=90, file=trim(dx_path))
    do i=1,ncomp-1
        read(90, *) x, dx(i)
    end do
    close(90)


    ! reading Strickler's coefficient at each section
    open(unit=85,file=trim(manning_strickler_path), status='unknown')
    do i=1,ncomp
        read(85, *) x, sk(i)
        call readXsection(i,(1.0/sk(i)),timesDepth)
        ! setting initial condition of water level
        oldY(i) = yy !+ z(i)
    end do
    close(85)


    !do i=1,ncomp
    !    call create_I2(i,ncomp)
    !end do


    ityp = 1


    ! setting initial condition of Q
    oldQ = qq

        ! setting initial condition
    ! setting initial condition from previous work
    !open(unit=91,file=trim(other_input)//'initialCondition.txt', status='unknown')
    !read(91, *)
    !do i=1,ncomp
    !    read(91, *) x, oldQ(i) !, oldY(i), oldArea(i)
    !end do
    !close(91)



    ! reading Q-Strickler's coefficient multiplier table
    do i=1,noQSKtable
        write(file_num,'(i4.4)')i
        open(86,file=trim(QSKtablePath)//'Q_Mannings_table_'//file_num//'.txt')
        do n=1,maxTableLength
            read(86,*,end=300) Q_Sk_Table(1, n, i), Q_Sk_Table(2, n, i)
        end do
300     close(86)
        Q_sk_tableEntry(i) = n-1
    end do

    x = 0.0


    ! Read hydrograph input Upstream
    open(unit=87, file=upstream_path)
    do n=1,boundaryFileMaxEntry
        read(87,*,end=301) USBoundary(1, n), USBoundary(2, n)
    end do
301 close(87)
    ppp = n-1   ! stores the number of entries in the boundary file


    ! Read hydrograph input Downstream
    open(unit=88, file=downstream_path)
    do n=1,boundaryFileMaxEntry
      read(88,*,end=302)  DSBoundary(1, n), DSBoundary(2, n)
    end do
302 close(88)
    qqq = n-1   ! stores the number of entries in the boundary file


   ! starting time
    !t = 28800.0
    !t = 315532800.0
    t = t0*60.0     !!! t0 is in hour. t is in minutes
    !t = 283996800.0 ! to starts from 2018
    !t = 1468800.0


    ! applying boundary
    ! interpolation of boundaries at the initial time step
    oldQ(1)    =r_interpol_time(USBoundary(1, 1:ppp),USBoundary(2, 1:ppp),ppp,t)
    oldY(ncomp)=r_interpol_time(DSBoundary(1, 1:qqq),DSBoundary(2, 1:qqq),qqq,t)


    ! DS Boundary treatment: from water level to area time series
    !ncompElevTable = xsec_tab(1,:,ncomp)
    !ncompAreaTable = xsec_tab(2,:,ncomp)

    xt=oldY(ncomp)
    !oldArea(ncomp)=r_interpol(ncompElevTable,ncompAreaTable,nel,xt)


    ! read lateral flow conditions
    do i=1,noLatFlow
        write(file_num,'(i4.4)')latFlowLocations(i)
        open(89,file=trim(lateralFlow_path)//'lateral_'//file_num//'.txt')
        do n=1,boundaryFileMaxEntry
            read(89,*,end=303) lateralFlowTable(1, n, i), lateralFlowTable(2, n, i)
        end do
303     close(89)
        dataInEachLatFlow(i) = n-1
    end do



    !bottomSlope = 9.23721E-05 !!! average slope of MS from BR to BC


    ! Converting initial WL to initial Area
    !do i = 1,ncomp-1
!       Now calculate y based on area calculated
!-------------------------------------
        !elevTable(:) = xsec_tab(1,:,i)
        !areaTable(:) = xsec_tab(2,:,i)

!       interpolate the cross section attributes based
        !oldArea(i)=r_interpol(elevTable,areaTable,nel,oldY(i))
!-------------------------------------
    !end do


    ! Open files for output
    path = trim(output_path) // 'output_wl.txt'
    open(unit=8, file=trim(path), status='unknown')

    path = trim(output_path) // 'q.txt'
    open(unit=9, file=trim(path), status='unknown')



    path = trim(output_path) // 'width.txt'
    open(unit=91, file=trim(path), status='unknown')

    path = trim(output_path) // 'area.txt'
    open(unit=92, file=trim(path), status='unknown')

    path = trim(output_path) // 'pere.txt'
    open(unit=93, file=trim(path), status='unknown')

    path = trim(output_path) // 'celerity.txt'
    open(unit=94, file=trim(path), status='unknown')

    !path = trim(output_path) // 'area.txt'
    !open(unit=51, file=trim(path), status='unknown')



    ! Output initial condition

    write(8, 10)  t*60., (oldY(i), i=1,ncomp)
    write(9, 10)  t*60., (oldQ(i), i=1, ncomp)
    !write(51, 10) t, (oldArea(i), i=1, ncomp)







    !bo = 10 ! Case specific condition

    !! Calculating initial values of diffusivity

    !celerity = 0.5
    !diffusivity = 10.


	theta = 1.0

    !pause 110
    !

    qpx = 0.
    !qp  = 0.


	!! some initial values
    width = 800. !!! average width of MS from BR to BC
    celerity = 1.0
    diffusivity = 3000.

    !do n=1, ntim-1
    do while ( t .lt. tfin *60.)

        ! new approach to add multiple lateral flow to the same node
        lateralFlow = 0
        do i=1,noLatFlow
            if (latFlowType(i) .eq. 1) then
                latFlowValue = r_interpol_time(lateralFlowTable(1, 1:dataInEachLatFlow(i), i), &
                    lateralFlowTable(2, 1:dataInEachLatFlow(i), i),dataInEachLatFlow(i),t)
            elseif (latFlowType(i) .eq. 2) then
                latFlowValue = r_interpol(lateralFlowTable(1, 1:dataInEachLatFlow(i), i), &
                    lateralFlowTable(2, 1:dataInEachLatFlow(i), i),dataInEachLatFlow(i),oldQ(latFlowLocations(i)-1))
            endif
            latFlowValue = latFlowValue / &
                    !sum(dx(latFlowLocations(i):latFlowLocations(i)+latFlowXsecs(i)-1))          !!! check this line
                    sum(dx(latFlowLocations(i)-1:latFlowLocations(i)-1+latFlowXsecs(i)-1))

            do j=1,latFlowXsecs(i)
                lateralFlow(latFlowLocations(i)+j-1)=lateralFlow(latFlowLocations(i)+j-1) + latFlowValue
            end do
        end do

        call calculateDT(t0, t,saveInterval, cfl, tfin)
        !print*, t, dtini
!!++++++++++++++++++++ Diffusive wave Forward sweep starts +++++++++++++++++++!!
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

            if (i .eq. ncomp) then
                alpha = 1.0
            else
                alpha = dx(i) / dx(i-1)
            end if


			qy   = a1 * oldQ(i-1) + a2 * oldQ(i) + a3 * qpx(i-1) + a4 * qpx(i)
			qxy  = b1 * oldQ(i-1) + b2 * oldQ(i) + b3 * qpx(i-1) + b4 * qpx(i)
			qxxy = dd1* oldQ(i-1) + dd2* oldQ(i) + dd3* qpx(i-1) + dd4* qpx(i)
			qxxxy= h1 * oldQ(i-1) + h2 * oldQ(i) + h3 * qpx(i-1) + h4 * qpx(i)


            ppi = - theta * diffusivity(i) * dtini / ( dx(i-1) ** 2.0 ) * 2.0 / (alpha*(alpha + 1.0)) * alpha
			qqi = 1.0 - ppi * (alpha + 1.0) / alpha
			rri = ppi / alpha
			ssi = qy  + dtini * diffusivity(i) * ( 1.0 - theta ) * qxxy + dtini * celerity(i) * lateralFlow(i)
			sxi = qxy + dtini * diffusivity(i) * ( 1.0 - theta ) * qxxxy

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
        qp(ncomp) = oldQ(ncomp-1)
        qpx(ncomp)= 0.

        do i = ncomp-1,1,-1

			qp(i) = eei(i) * qp(i+1) + ffi(i)
			qpx(i)= exi(i) *qpx(i+1) + fxi(i)
        end do


        qp(1) =r_interpol_time(USBoundary(1, 1:ppp),USBoundary(2, 1:ppp),ppp,t+dtini/60.)
        newY(ncomp) =r_interpol_time(DSBoundary(1, 1:qqq),DSBoundary(2, 1:qqq),qqq,t+dtini/60.)
        depth(ncomp)=newY(ncomp)-z(ncomp)

        do i=ncomp,1,-1

            q_sk_multi = 1.0
            do pp = 1, size(Q_sk_tableEntry)
                if (  ( eachQSKtableNodeRange(1,pp) - i) * ( eachQSKtableNodeRange(2,pp) - i) .le. 0 ) then
                    q_sk_multi = r_interpo_nn(Q_Sk_Table(1,1:Q_sk_tableEntry(pp),pp),   &
                      Q_Sk_Table(2,1:Q_sk_tableEntry(pp),pp),Q_sk_tableEntry(pp),qp(i))
                end if
            end do

    !      Nazmul change: read all attributes from tab file
            elevTable = xsec_tab(1,:,i)
            convTable = xsec_tab(5,:,i)

            areaTable = xsec_tab(2,:,i)
            pereTable = xsec_tab(3,:,i)
            topwTable = xsec_tab(6,:,i)
    !     interpolate the cross section attributes based on water elevation
            xt=newY(i)
            co(i)  =q_sk_multi * r_interpol(elevTable,convTable,nel,xt)
            width = r_interpol(elevTable,topwTable,nel,xt)



            area(i) = r_interpol(elevTable,areaTable,nel,xt)
            pere(i) = r_interpol(elevTable,pereTable,nel,xt)
            bo(i) = r_interpol(elevTable,topwTable,nel,xt)


            sfi = ( qp(i) / co(i) ) ** 2.0

            celerity2(i)=5.0 / 3.0 * sfi ** 0.3 * abs(qp(i)) ** 0.4 / width ** 0.4 / (1/(sk(i)*q_sk_multi)) ** 0.6
            diffusivity2(i) = abs(qp(i)) / 2.0 / width / sfi


            if (i .gt. 1) newY(i-1) = newY(i) + sign ( sfi, qp(i) ) * dx(i-1)
        end do

        do i=1,ncomp
            !celerity(i) = sign ( sum(celerity2) / ncomp, qp(i) )

            celerity(i) =  sum(celerity2) / ncomp
            !print*, celerity(i)
        end do


        diffusivity = sum(diffusivity2) / ncomp

		do i = 1, ncomp
			if (diffusivity(1) .gt. 10) diffusivity = 10
		end do
        ! Final update

        newQ = qp


        !diffusivity = 0.0


        !celerity = 0.5
        !diffusivity = 10.

        !pause 120


        !! Saving results




        t = t + dtini / 60.
        !print "('- cycle',i9,'  completed')", n
        print "('- simulation time ',f16.2,'  completed')", t




        !! Saving results
        if ( (mod( t-t0*60.  ,saveInterval/60. ) .eq. 0.0) .or. ( t .eq. tfin *60. ) ) then
            write(8, 10) t*60., (newY(i), i=1,ncomp)
            write(9, 10) t*60., (newQ(i), i=1,ncomp)


          !  write(91, 10) t*60., (bo(i), i=1,ncomp)
          !  write(92, 10) t*60., (area(i), i=1,ncomp)
          !  write(93, 10) t*60., (pere(i), i=1,ncomp)
          !  write(94, 10) t*60., (celerity(i), i=1,ncomp)


            !write(51, 10) t*60., (lateralFlow(i), i=1, ncomp)
            !write(51, 10) t*60., (diffusivity2(i), i=1, ncomp)
        end if

        oldQ   = newQ
        oldY = newY
        !print*, newQ
    !pause 120
    end do
    ! End of time loop


    close(8)
    close(9)
    close(51)
    close(91)
    close(92)
    close(93)
    close(94)
    !close(81)

    print*, 'dx', (dx(i), i=1, ncomp-1)
    !print*, 'Froude', (froud(i), i=1, ncomp)
    print*, 'Bed', (z(i), i=1, ncomp)
    !print*, 'newArea', (newArea(i), i=1, ncomp)
    !print*, 'Courant no', (courant(i), i=1, ncomp-1)
    !print*, 'Maximum Courant no', maxCourant

    !
!10  format(f12.2 , <ncomp>f13.3)
!11  format(f12.2, <ncomp>f13.9)
10  format(f12.2 , 1500f13.3)
11  format(f12.2, 1500f13.9)

    call cpu_time( t2 )
    print '("Time = ",f10.3," seconds.")',t2 - t1
    !pause 202
end program mesh


function r_interpol(x,y,jj,xt)

    integer, intent(in) :: jj
    real, intent(in) :: xt, x(jj), y(jj)

    if (xt.lt.maxval(x) .and. xt.ge.minval(x)) then
        do j=1,jj-1
            if((x(j)-xt)*(x(j+1)-xt).le.0)then

                yt=(xt-x(j))/(x(j+1)-x(j))*(y(j+1)-y(j))+y(j)

                EXIT
            endif
        end do
    else
        print*, xt, ' is not within the limit'
        print*, 'maxval(x)= ', maxval(x), 'and minval(x)=', minval(x),'so',  xt, ' is not within the limit'
        print*, 'jj', jj
        print*, 'x', (x(i), i=1, jj)
        print*, 'y', (y(i), i=1, jj)
        stop
        !if (xt.le. minval(x)) yt=minval(y)
        !if (xt.ge. maxval(x)) yt=maxval(y)
    end if
    r_interpol = yt
    return
end function

function r_interpol_time(x,y,jj,xt)

    integer(kind=4), intent(in) :: jj
    real(kind=4), intent(in) :: x(jj), y(jj)
    real(kind=4), intent(in) :: xt
    real(kind=4) :: yt
    !real(kind=8), intent(out) :: r_interpol_time


    if (xt.lt.maxval(x) .and. xt.ge.minval(x)) then
        do j=1,jj-1
            if((x(j)-xt)*(x(j+1)-xt).le.0)then

                yt=(xt-x(j))/(x(j+1)-x(j))*(y(j+1)-y(j))+y(j)

                EXIT
            endif
        end do
    else
        print*, xt, ' is not within the limit'
        print*, 'maxval(x)= ', maxval(x), 'and minval(x)=', minval(x),'so',  xt, ' is not within the limit'
        print*, 'jj', jj
        print*, 'x', (x(i), i=1, jj)
        print*, 'y', (y(i), i=1, jj)
        stop
        !if (xt.le. minval(x)) yt=minval(y)
        !if (xt.ge. maxval(x)) yt=maxval(y)
    end if
    r_interpol_time = yt
    ! print*,xt
    return
end function


function r_interpo_nn(x,y,jj,xt)

    integer, intent(in) :: jj
    real, intent(in) :: xt, x(jj), y(jj)

    real :: yt

    ! nn means nearest neighbour

    if (xt.le. x(1)) then
        yt=y(1)
    elseif (xt.ge. x(jj)) then
        yt=y(jj)
    else
        do j=1,jj-1
            if((x(j)-xt)*(x(j+1)-xt).le.0)then

                yt=(xt-x(j))/(x(j+1)-x(j))*(y(j+1)-y(j))+y(j)

                EXIT
            endif
        end do
    end if
    r_interpo_nn = yt
    return
end function
