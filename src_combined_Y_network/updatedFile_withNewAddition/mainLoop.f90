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
    integer(kind=4) :: i, j, k, ppn, qqn, n, ntim, igate, pp, boundaryFileMaxEntry, saveFrequency, maxNotSwitchRouting


    real(kind=4) :: cour, da, dq, yy, x, saveInterval, width
    real(kind=4) :: qq, qn, xt, r_interpol, maxCourant, dtini_given
    real(kind=4) :: arean, areac, hyrdn, hyrdc, perimn, perimc, qcrit, s0ds, timesDepth, latFlowValue

    real(kind=4) :: t, r_interpol_time, tfin, t1, t2, t0 !t0 start time
    !doubleprecision ::  t, r_interpol_time, tfin, t1, t2, t0 !t0 start time

    integer(kind=4) :: tableLength
    real(kind=4) :: area_0, width_0, errorY, hydR_0, q_sk_multi

    real(kind=4) :: r_interpo_nn

    character(len=128) :: output_path, other_input
    character(len=128) :: path

    !real(kind=4) :: pere(500)


    ! open file for input data

    call cpu_time( t1 )

    !open(unit=1,file="../lower_Mississippi/input/input_BR_2_BC_2009.txt",status='unknown')
    !open(unit=1,file="../lower_Mississippi/input/input_BR_2_BC_2009.txt",status='unknown')
    !open(unit=1,file="../lateralFlow_test/input/input_dynamic_lateralFlow.txt",status='unknown')
    !open(unit=1,file="../Mississippi_River_11_years_20200511/input/input_Mississippi_BR2SWP_dynamic_new.txt",status='unknown')
    !open(unit=1,file="../Vermelion_River/input/input_Vermelion_dynamic_20200526.txt",status='unknown')
    !open(unit=1,file="../../../MESH_code/4-A-1_US_2bound/input_naturalChannel_tidal.txt",status='unknown')
    open(unit=1,file="../Rectangular_Y_Channel/input/input_naturalChannel.txt",status='unknown')

    print*, 'Reading input file'

    ! read data
    read(1,*) dtini_given     ! in seconds
    dtini = dtini_given         !; print*, dtini; pause 500
    read(1,*) dxini
    read(1,*) t0        ! in hours
    read(1,*) tfin      ! in hours
	ntim = floor( (tfin - t0) / dtini * 3600)
	print*, ntim
	read(1,*) nlinks
	allocate(nx1(nlinks))
	read(1,*) (nx1(i), i=1, nlinks)
    read(1,*) phi
    read(1,*) theta
    read(1,*) thetas
    read(1,*) thesinv
    read(1,*) alfa2
    read(1,*) alfa4
    read(1,*) f
    read(1,*) skk
    read(1,*) yy !; print*, yy; pause 5000
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


    allocate(notSwitchRouting(nlinks))
    allocate(currentROutingDiffusive(nlinks))
    allocate(xSection_path(nlinks))
    do i=1,nlinks
        read(1,*) xSection_path(i)
    end do

    allocate(manning_strickler_path(nlinks))
    do i=1,nlinks
        read(1,*) manning_strickler_path(i)
    end do

    read(1,*) nupbnds
    allocate(upBoundTableEntry(nupbnds))
    allocate(upstream_path(nupbnds))
    do i=1,nupbnds
        read(1,*) upstream_path(i)
    end do

    read(1,*) ndnbnds
    allocate(downBoundTableEntry(ndnbnds))
    allocate(downstream_path(ndnbnds))
    do i=1,ndnbnds
        read(1,*) downstream_path(i)
    end do

    allocate(QSKtablePath(nlinks))
    do i=1,nlinks
        read(1,*) QSKtablePath(i)
    end do

    allocate(dx_path(nlinks))
    do i=1,nlinks
        read(1,*) dx_path(i)
    end do

    allocate(lateralFlow_path(nlinks))
    do i=1,nlinks
        read(1,*) lateralFlow_path(i)
    end do
    read(1,*) output_path
    read(1,*) option_dsbc
    read(1,*) maxTableLength
    read(1,*) nel
    read(1,*) timesDepth
    read(1,*) other_input
    read(1,*) boundaryFileMaxEntry
    read(1,*) saveInterval !; saveFrequency = saveInterval / dtini

    allocate(noLatFlow(nlinks))
	read(1,*) (noLatFlow(i), i=1, nlinks)

    allocate(latFlowLocations(maxval(noLatFlow),nlinks))   ! all the first nodes where a lateral flow starts
    allocate(latFlowType(maxval(noLatFlow),nlinks))        ! Lateral flow type: Type = 1 = time series; Type 2 = flow as a function of upstream flow
    allocate(latFlowXsecs(maxval(noLatFlow),nlinks))       ! no of x-secs at the downstream that the lateral flow is applied

    do j = 1,nlinks
        ncomp=nx1(j)
        read(1,*) (latFlowLocations(i,j), i=1, noLatFlow(j))
        do i=1,noLatFlow(j)
            if ((latFlowLocations(i,j)-1)*(latFlowLocations(i,j)-ncomp) .eq. 0) then
                print*, 'ERROR: Lateral flow cannot be applied at the boundaries'
                stop
            end if
        end do
    end do
    do j = 1,nlinks
        ncomp=nx1(j)
        read(1,*) (latFlowType(i,j), i=1, noLatFlow(j))
        do i=1,noLatFlow(j)
            if (latFlowType(i,j) .eq. 1) then
                print*, 'Lateral flow at node = ', latFlowLocations(i,j), ', is a time series at reach ', j
            elseif (latFlowType(i,j) .eq. 2) then
                print*, 'Lateral flow at node = ', latFlowLocations(i,j), ', is a function of upstream flow at reach ', j
            else
                print*, 'Wrong lateral flow type is provided. Type ', latFlowType(i,j), 'is not a valid type at reach ', j
                stop
            end if
        end do
    end do

    do j = 1,nlinks
        ncomp = nx1(j)
        read(1,*) (latFlowXsecs(i,j), i=1, noLatFlow(j))
    end do


	allocate(noQSKtable(nlinks))                                ! how many tables are there in each river reach
	read(1,*) (noQSKtable(i), i=1, nlinks)

    allocate(eachQSKtableNodeRange(2,maxval(noQSKtable),nlinks))
    do j = 1,nlinks
        read(1,*) (eachQSKtableNodeRange(1,i,j), i=1, noQSKtable(j))    ! upper limit of the river node number assigned to current table
        read(1,*) (eachQSKtableNodeRange(2,i,j), i=1, noQSKtable(j))    ! lower limit of the river node number assigned to current table
    !!! Need to test so that one section does not corresponds to more than one table
        do i = 2, noQSKtable(j)
            if ( eachQSKtableNodeRange(2,i-1,j) .ge. eachQSKtableNodeRange(1,i,j) ) then
                print*, 'Wrong range of nodes applied for Q-Sk table.'
                print*, 'Lower limit of Table ', i-1,'must be smaller than the upper limit of Table ', i, ' of reach', j
                stop
            end if
        end do
    end do
    close (1)

    ! Allocate arrays
    call setup_arrays(ntim, maxval(nx1), maxTableLength, boundaryFileMaxEntry, maxval(noLatFlow), maxval(noQSKtable), nlinks)
    call setup_arrays_section
    call setup_xsec_attribute_module(nel, maxval(nx1),nlinks)

    dt = dtini

    do j = 1,nlinks
        ncomp = nx1(j)
        open(unit=90, file=trim(dx_path(j)), status='unknown')
            do i=1,ncomp-1
                read(90, *) x, dx(i,j)
            end do
        close(90)
    end do

    ! reading Strickler's coefficient at each section
    do j = 1,nlinks
        ncomp = nx1(j)
        open(unit=85,file=trim(manning_strickler_path(j)), status='unknown') !! //'Mannings_Stricklers_coeff.txt', status='unknown')
        do i=1,ncomp
            read(85, *) x, sk(i,j)
            call readXsection(i,(1.0/sk(i,j)),timesDepth,j)
            ! This subroutine creates attribute table for each cross sections and saves in the hdd
            ! setting initial condition
            oldY(i,j) = yy  !+ z(i,j)
        end do
        close(85)
    end do

    do j = 1,nlinks
        ncomp = nx1(j)
        do i=1,ncomp
            call create_I2(i,ncomp,j)
        end do
    end do
    NAnum = -100
    ityp = 1

    ! setting initial condition
    ! setting initial condition from previous work
    !open(unit=91,file=trim(output_path)//'initialCondition.txt', status='unknown')
    ! read(91, *)
    !do i=1,ncomp
    !    read(91, *) oldQ(i), oldY(i), oldArea(i)
    !end do
    !close(91)

    !q(1, :) = qq
    oldQ = qq

    ! reading Q-Strickler's coefficient multiplier table
    do j = 1,nlinks
        ncomp = nx1(j)
        do i=1,noQSKtable(j)
            write(file_num,'(i4.4)')i
            open(86,file=trim(QSKtablePath(j))//'Q_Mannings_table_'//file_num//'.txt')
            do n=1,maxTableLength
                read(86,*,end=300) Q_sk_Table(1, n, i, j), Q_sk_Table(2, n, i, j)
            end do
300         close(86)
            Q_sk_tableEntry(i,j) = n-1
        end do
    end do



    !!! Code from DongHa

    !* the number of links that are immediately upstream of link j
    ndep(1)=0; ndep(2)=0; ndep(3)=2

    !* link number of k_th link that is immediately upstream of link j
    uslinks(1,1)=NAnum  !not available
    uslinks(1,2)=NAnum  !not available
    uslinks(1,3)=1;  uslinks(2,3)=2
    !* link number that is immediately downstream of link j
    dslink(1)=3; dslink(2)=3; dslink(3)=NAnum

    !* when data at either upper or lower end of link j is available,
    !* instrdflag(j,1)=1 when water level data is known at the upper end of link j
    !* instrdflag(j,1)=2 when discharge data is known
    !* instrdflag(j,1)=3 when rating curve is known

    !* instrdflag(j,2)=1 when water level data is known at the lower end of link j
    !* instrdflag(j,2)=2 when discharge data is known
    !* instrdflag(j,2)=3 when rating curve is known
    !* Otherwise, instrdflag(j,1/2)=0
    !! Nazmul: This part is currently hard coded, but we need to move it in a flexible manner.

    print*, 'Nazmul 1'
    instrdflag(1,1)=2; instrdflag(1,2)=1    !*discharge data is known at the upper end of link 1
    instrdflag(2,1)=2; instrdflag(2,2)=1    !*discharge is known at the upper end of link 2
    instrdflag(3,1)=2; instrdflag(3,2)=1    !*stage data is known at the lower end of link 3

    ! Nazmul: Need to create a connectivity table like the following:
    ! p=Junction sequence no, q= number of river reaches connected to that junction, (riverReachSequenceNoInThatJunction(i),i=1,q), (streamOrderOfEachRivers(i),i=1,q)



    print*, 'Nazmul, 11'
    x = 0.0

    ! Read hydrograph input Upstream
    do j = 1, nupbnds
        open(unit=87, file=upstream_path(j))
        do n=1,boundaryFileMaxEntry
            read(87,*,end=301) USBoundary(1, n,j), USBoundary(2, n,j)       !! time column is in minutes
        end do
301     close(87)
        upBoundTableEntry(j) = n-1   ! stores the number of entries in the boundary file
    end do

    ! Read hydrograph input Downstream
    do j = 1, ndnbnds
        open(unit=88, file=downstream_path(j))
        do n=1,boundaryFileMaxEntry
          read(88,*,end=302)  DSBoundary(1, n, j), DSBoundary(2, n, j)       !! time column is in minutes
        end do
302     close(88)
        downBoundTableEntry(j) = n-1   ! stores the number of entries in the boundary file
    end do

    t=t0*60.0     !! from now on, t is in minute










    ! applying boundary
    ! interpolation of boundaries at the initial time step
    !! Need to define which channel has terminal boundary
    do j = 1, nlinks
        ncomp = nx1(j)
        ppn = 1
        if (instrdflag(j,1) .eq. 2) then
            ! interpolation of boundaries at the desired time step
            oldQ(1,j)=r_interpol_time(USBoundary(1, 1:upBoundTableEntry(ppn), ppn), &
                USBoundary(2, 1:upBoundTableEntry(ppn), ppn),upBoundTableEntry(ppn),t)
            ppn = ppn +1
        end if
        !print*, t+dtini, newQ(1)
        qqn = 1
        if (instrdflag(j,2) .eq. 1) then
            ! interpolation of boundaries at the desired time step
            oldY(ncomp,j)=r_interpol_time(DSBoundary(1, 1:downBoundTableEntry(qqn), qqn), &
                DSBoundary(2, 1:downBoundTableEntry(qqn), qqn),downBoundTableEntry(qqn),t)
            qqn = qqn +1
            ncompElevTable = xsec_tab(1,:,ncomp,j)
            ncompAreaTable = xsec_tab(2,:,ncomp,j)
            xt=oldY(ncomp,j)
            oldArea(ncomp,j)=r_interpol(ncompElevTable,ncompAreaTable,nel,xt)


        end if
        !print*, j, ncomp, instrdflag(j,2), oldY(ncomp,j), oldArea(ncomp,j)
    end do


    ! Applying initial condition of area
    do j=1, nlinks
        ncomp = nx1(j)
        do i=1,ncomp-1
            elevTable = xsec_tab(1,1:nel,i,j)
            areaTable = xsec_tab(2,1:nel,i,j)
            oldArea(i,j)=r_interpol(elevTable,areaTable,nel,oldY(i,j))
        enddo
    end do


    ! read lateral flow conditions
    do j=1, nlinks
        ncomp = nx1(j)
        do i=1,noLatFlow(j)
            write(file_num,'(i4.4)')latFlowLocations(i,j)
            open(89,file=trim(lateralFlow_path(j))//'lateral_'//file_num//'.txt')
            do n=1,boundaryFileMaxEntry
                read(89,*,end=303) lateralFlowTable(1, n, i,j), lateralFlowTable(2, n, i,j)
            end do
303         close(89)
            dataInEachLatFlow(i,j) = n-1
        end do
    end do


    ! Open files for output
    path = trim(output_path) // 'output_wl.txt'
    open(unit=8, file=trim(path), status='unknown')
    path = trim(output_path) // 'q.txt'
    open(unit=9, file=trim(path), status='unknown')

    path = trim(output_path) // 'area.txt'
    open(unit=51, file=trim(path), status='unknown')

    path = trim(output_path) // 'width.txt'
    open(unit=91, file=trim(path), status='unknown')

    path = trim(output_path) // 'pere.txt'
    open(unit=93, file=trim(path), status='unknown')

    path = trim(output_path) // 'dimensionless_Cr.txt'
    open(unit=941, file=trim(path), status='unknown')

    path = trim(output_path) // 'dimensionless_Fr.txt'
    open(unit=942, file=trim(path), status='unknown')

    path = trim(output_path) // 'dimensionless_Fi.txt'
    open(unit=943, file=trim(path), status='unknown')

    path = trim(output_path) // 'dimensionless_Di.txt'
    open(unit=944, file=trim(path), status='unknown')

    path = trim(output_path) // 'dimensionless_Fc.txt'
    open(unit=945, file=trim(path), status='unknown')

    path = trim(output_path) // 'dimensionless_D.txt'
    open(unit=946, file=trim(path), status='unknown')

    path = trim(output_path) // 'celerity.txt'
    open(unit=95, file=trim(path), status='unknown')

    path = trim(output_path) // 'diffusivity.txt'
    open(unit=96, file=trim(path), status='unknown')

    path = trim(output_path) // 'routing.txt'
    open(unit=97, file=trim(path), status='unknown')


	! Some essential initial parameters for Diffusive Wave
	theta = 1.0
	qpx = 0.

    width = 10. !!! average width of MS from BR to BC
    celerity = 1.0
    diffusivity = 10.


    !!! setting initial values of dimensionless parameters
    !dimensionless_Cr, dimensionless_Fo, dimensionless_Fi, dimensionless_Fc, dimensionless_Di, dimensionless_D
    dimensionless_Fi = 10.1
    dimensionless_Fc = 10.1
    currentROutingDiffusive = 0





        ! Output initial condition
    do j=1, nlinks
        ncomp = nx1(j)

        write(8, 10)  t*60.0, j, (oldY(i,j), i=1,maxval(nx1))
        write(9, 10)  t*60.0, j, (oldQ(i,j), i=1, ncomp)
        write(51, 10) t*60.0, j, (oldArea(i,j), i=1, ncomp)
        !write(91, 10) t*60.0, j, (bo(i), i=1, ncomp)
        !write(93, 10) t*60.0, j, (pere(i), i=1, ncomp)
        write(941, 10) t*60.0, j, (dimensionless_Cr(i), i=1, ncomp-1)
        write(942, 10) t*60.0, j, (dimensionless_Fo(i), i=1, ncomp-1)
        write(943, 10) t*60.0, j, (dimensionless_Fi(i), i=1, ncomp-1)
        write(944, 10) t*60.0, j, (dimensionless_Di(i), i=1, ncomp-1)
        write(945, 10) t*60.0, j, (dimensionless_Fc(i), i=1, ncomp-1)
        write(946, 10) t*60.0, j, (dimensionless_D(i), i=1, ncomp-1)
        write(95, 10) t*60.0, j, (celerity2(i), i=1, ncomp)
        write(96, 10) t*60.0, j, (diffusivity2(i), i=1, ncomp)
        write(97, *) t*60.0, j, currentROutingDiffusive(j)

    end do


    frus2 = 9999.
    notSwitchRouting=0
    maxNotSwitchRouting = 10
    !
    ! Loop on time
    !
    !do n=0, ntim-1
    do while ( t .lt. tfin *60.)
    !+++---------------------------------------------------------------------------------+
    !+ Run predictor step though all the links
    !+++---------------------------------------------------------------------------------+
    ppn = 1
    qqn = 1

    do j = 1, nlinks

        ncomp = nx1(j)
        !+++---------------------------------------------------------------------------------+
        !+ Run predictor step though all the links
        !+++---------------------------------------------------------------------------------+

        if (instrdflag(j,1) .eq. 2) then
            ! interpolation of boundaries at the desired time step
            newQ(1,j)=r_interpol_time(USBoundary(1, 1:upBoundTableEntry(ppn), ppn), &
                USBoundary(2, 1:upBoundTableEntry(ppn), ppn),upBoundTableEntry(ppn),t+dtini/60.)
            ppn = ppn +1
        end if
        !print*, t+dtini, newQ(1)

        if (instrdflag(j,2) .eq. 1) then
            ! interpolation of boundaries at the desired time step
            newY(ncomp,j)=r_interpol_time(DSBoundary(1, 1:downBoundTableEntry(qqn), qqn), &
                DSBoundary(2, 1:downBoundTableEntry(qqn), qqn),downBoundTableEntry(qqn),t+dtini/60.)

            ncompElevTable = xsec_tab(1,:,ncomp,j)
            ncompAreaTable = xsec_tab(2,:,ncomp,j)
            xt=newY(ncomp,j)
            newArea(ncomp,j)=r_interpol(ncompElevTable,ncompAreaTable,nel,xt)

           ! print*, 'downBoundTableEntry(qqn)', qqn, downBoundTableEntry(qqn), DSBoundary(1, 1:downBoundTableEntry(qqn), qqn), &
            !    DSBoundary(2, 1:downBoundTableEntry(qqn), qqn), j, newY(ncomp,j)
            qqn = qqn +1
        end if







          !! Calculating Y normal at the upstream!!
        S_0 = (-z(2,j)+z(1,j))/dx(1,j)
        if (S_0 .gt. 0.) then
            q_sk_multi = 1.0
            do pp = 1, noQSKtable(j)
                if (  ( eachQSKtableNodeRange(1,pp,j) - 1) * ( eachQSKtableNodeRange(2,pp,j) - 1) .le. 0 ) then
                    tableLength = Q_sk_tableEntry(pp,j)
                    q_sk_multi = r_interpo_nn(Q_sk_Table(1,1:tableLength,pp,j), &
                        Q_sk_Table(2,1:tableLength,pp,j),tableLength,newQ(1,j))
                end if
            end do

            elevTable = xsec_tab(1,:,1,j)
            areaTable = xsec_tab(2,:,1,j)
            rediTable = xsec_tab(4,:,1,j)
            topwTable = xsec_tab(6,:,1,j)
            area_0 = r_interpol(elevTable,areaTable,nel,oldY(1,j))
            width_0= r_interpol(elevTable,topwTable,nel,oldY(1,j))
            area_crit=area_0
            errorY = 100.
            do while (errorY .gt. 0.0001)

                hydR_0 = r_interpol(areaTable,rediTable,nel,area_0)
                area_norm = newQ(1,j)/sk(1,j)/q_sk_multi/ hydR_0 ** (2./3.) / sqrt(S_0)

                errorY = abs(area_norm - area_0) / area_0
                area_0 = area_norm
                area_crit= (newQ(1,j) * newQ(1,j) * width_0 / grav) ** (1./3.)
                width_0  = r_interpol(areaTable,topwTable,nel,area_crit)

            enddo
            y_norm_us = r_interpol(areaTable,elevTable,nel,area_0)
            y_crit_us = r_interpol(areaTable,elevTable,nel,area_crit)
        endif

        ! new approach to add multiple lateral flow to the same node
        lateralFlow = 0
        do i=1,noLatFlow(j)
            if (latFlowType(i,j) .eq. 1) then
                latFlowValue = r_interpol_time(lateralFlowTable(1, 1:dataInEachLatFlow(i,j), i,j), &
                    lateralFlowTable(2, 1:dataInEachLatFlow(i,j), i,j),dataInEachLatFlow(i,j),t)
                !print*, t, latFlowValue
            elseif (latFlowType(i,j) .eq. 2) then
                latFlowValue = r_interpol(lateralFlowTable(1, 1:dataInEachLatFlow(i,j), i,j), &
                    lateralFlowTable(2, 1:dataInEachLatFlow(i,j), i,j),dataInEachLatFlow(i,j),oldQ(latFlowLocations(i,j)-1,j))
            endif
            latFlowValue = latFlowValue / &
                    !sum(dx(latFlowLocations(i):latFlowLocations(i)+latFlowXsecs(i)-1))          !!! check this line
                    sum(dx(latFlowLocations(i,j)-1:latFlowLocations(i,j)-1+latFlowXsecs(i,j)-1,j))

            do k=1,latFlowXsecs(i,j)
                lateralFlow(latFlowLocations(i,j)+k-1)=lateralFlow(latFlowLocations(i,j)+k-1) + latFlowValue
            end do

        end do

        lowerLimitCount = 0; higherLimitCount = 0

        do i=1,ncomp-1
            if ((dimensionless_Fi(i) .ge. 5.) .or. (dimensionless_Fc(i)  .ge. 5.))then
                higherLimitCount(j) = higherLimitCount(j) + 1
            elseif ((dimensionless_Fi(i) .le. 3.) .or. (dimensionless_Fc(i)  .le. 3.))then
                lowerLimitCount(j) = lowerLimitCount(j) + 1
            end if
        end do

        higherLimitCount = 0; lowerLimitCount = ncomp
            !print*, lowerLimitCount, higherLimitCount
            !if ( (maxval(dimensionless_Fi) .ge. 10.) .or.  &
            !     (maxval(dimensionless_Fc) .ge. 10.))then

            ! old switching algorithm
            !if (higherLimitCount .ge. ncomp/2.) then
            !    call mesh_diffusive(ppp,qqq, t0, t, tfin, saveInterval)
            !    currentROutingDiffusive = 1
            !elseif (lowerLimitCount .ge. ncomp/2.) then
            !    call mesh_dynamic(dtini_given, ppp,qqq, t0, t, tfin, saveInterval)
            !    currentROutingDiffusive = 0
            !else
            !    if (currentROutingDiffusive .eq. 1) then
            !        call mesh_diffusive(ppp,qqq, t0, t, tfin, saveInterval)
            !        currentROutingDiffusive = 1
            !    else
            !        call mesh_dynamic(dtini_given, ppp,qqq, t0, t, tfin, saveInterval)
            !        currentROutingDiffusive = 0
            !    end if
            !end if



            ! new switching algorithm
            ! not used for now
            ! for now the code will run dynamic routing for all the river reaches
        higherLimitCount(j) = ncomp
        if (higherLimitCount(j) .ge. ncomp/2.) then
            if ( (currentROutingDiffusive(j) .eq. 0) .and. (notSwitchRouting(j) .le. maxNotSwitchRouting)) then
                call mesh_dynamic_predictor(dtini_given, t0, t, tfin, saveInterval,j)
                currentROutingDiffusive(j) = 0
            else
    !            call mesh_diffusive(ppp,qqq, t0, t, tfin, saveInterval)
                if (currentROutingDiffusive(j) .eq. 0) notSwitchRouting(j) = 0
                currentROutingDiffusive(j) = 1
            end if
        elseif (lowerLimitCount(j) .ge. ncomp/2.) then
            if ( (currentROutingDiffusive(j) .eq. 1) .and. (notSwitchRouting(j) .le. maxNotSwitchRouting)) then
     !           call mesh_diffusive(ppp,qqq, t0, t, tfin, saveInterval)
                currentROutingDiffusive(j) = 1
            else
                call mesh_dynamic_predictor(dtini_given, t0, t, tfin, saveInterval,j)
                if (currentROutingDiffusive(j) .eq. 1) notSwitchRouting(j) = 0
                currentROutingDiffusive(j) = 0
            end if
        else
            if (currentROutingDiffusive(j) .eq. 1) then
     !           call mesh_diffusive(ppp,qqq, t0, t, tfin, saveInterval)
                currentROutingDiffusive(j) = 1
            else
                call mesh_dynamic_predictor(dtini_given, t0, t, tfin, saveInterval,j)
                currentROutingDiffusive(j) = 0
            end if
        end if

        !print*, newArea(35,1),'6'
    end do  ! end off j loop

    do j =  nlinks,-1,1
        ncomp = nx1(j)
        if (currentROutingDiffusive(j) .eq. 0) then
            call mesh_dynamic_corrector(dtini_given, t0, t, tfin, saveInterval,j)
        end if
    end do

        !print*, newArea(35,1),'7'
    !!! Calculate Fc and Fi
    !call calc_dimensionless_numbers

    do j = 1, nlinks
        ncomp = nx1(j)
        do i=1,ncomp
            froud(i)=abs(newQ(i,j))/sqrt(grav*newArea(i,j)**3.0/bo(i))
            if (i .lt. ncomp) then
                courant(i)=(newQ(i,j)+newQ(i+1,j))/(newArea(i,j)+newArea(i+1,j))*dtini/dx(i,j)
            endif
        enddo
    enddo

    if (maxCourant .lt. maxval (courant)) then
        maxCourant = maxval (courant)
    endif

    t = t + dtini/60.
    notSwitchRouting(j) = notSwitchRouting(j) + 1
    !t = t0*60. + (n+1)*dtini/60.
    !print "('- cycle',i9,'  completed')", n
    !print "('- simulation time ',2f16.2,'  completed')", t


    !if(mod(n+1,24*saveFrequency) .eq. 0 .or. (n.eq.0))write(*,*)'Nstep =', n, 'Days = ', t/60./24., real(n+1)/real(ntim)*100., '% completed'
    !print*, 'dqp', (dqp(i), i=1, ncomp)
    !print*, 'dqc', (dqc(i), i=1, ncomp)
    !print*, 'qp', (qp(i), i=1, ncomp)
    !print*, 'oldQ', (oldQ(i), i=1, ncomp)
    !print*, 'newQ', (newQ(i), i=1, ncomp) celerity

    !if (mod(n+1,saveFrequency) .eq. 0 .or. n .eq. (ntim-1)) then

    print*, 'times', t, t0, dtini, (currentROutingDiffusive(j),j=1,nlinks)
    !if ( (mod( (t-t0*60.)*60.  ,saveInterval) .eq. 0.0) .or. ( t .eq. tfin *60. ) ) then

    do j = 1, nlinks
        ncomp = nx1(j)
        write(8, 10) t*60.,j, (newY(i,j), i=1,41)
        write(9, 10) t*60.,j, (newQ(i,j), i=1,ncomp)
        write(51, 10) t*60.,j, (newArea(i,j), i=1, ncomp)

        write(91, 10) t*60.,j, (bo(i), i=1, ncomp)
        write(93, 10) t*60.,j, (pere(i), i=1, ncomp)

        write(941, 10) t*60.0,j, (dimensionless_Cr(i), i=1, ncomp-1)
        write(942, 10) t*60.0,j, (dimensionless_Fo(i), i=1, ncomp-1)
        write(943, 10) t*60.0,j, (dimensionless_Fi(i), i=1, ncomp-1)
        write(944, 10) t*60.0,j, (dimensionless_Di(i), i=1, ncomp-1)
        write(945, 10) t*60.0,j, (dimensionless_Fc(i), i=1, ncomp-1)
        write(946, 10) t*60.0,j, (dimensionless_D(i),  i=1, ncomp-1)

        !write(941, 10) t*60.0, (dimensionless_Cr(i), i=1, numScallingParameters)
        !write(942, 10) t*60.0, (dimensionless_Fo(i), i=1, numScallingParameters)
        !write(943, 10) t*60.0, (dimensionless_Fi(i), i=1, numScallingParameters)
       ! write(944, 10) t*60.0, (dimensionless_Di(i), i=1, numScallingParameters)
       ! write(945, 10) t*60.0, (dimensionless_Fc(i), i=1, numScallingParameters)
        !write(946, 10) t*60.0, (dimensionless_D(i), i=1, numScallingParameters)

        write(95, 10) t*60.,j, (celerity2(i), i=1, ncomp)
        write(96, 10) t*60.,j, (diffusivity2(i), i=1, ncomp)
        write(97, *) t*60.0,j, currentROutingDiffusive(j), notSwitchRouting(j)
    end do

    !end if
    !  if (n .eq. 8000) pause

    !write(81, *) t, newArea(ncomp)

    ! update of Y, Q and Area vectors
    oldY   = newY
    newY=0.0
    oldQ   = newQ
    newQ=0.0
    oldArea= newArea
    newArea=0.0

    !!! test
    !if (t .gt. 50000) dtini = 103
    !!! test end
    !pause 5000
    end do  ! end of time loop
    ! End of time loop

    close(8)
    close(9)
    close(51)

    close(91)
    close(93)
    close(941)
    close(942)
    close(943)
    close(944)
    close(945)
    close(946)
    close(95)
    close(96)
    !close(81)

    print*, 'dx', (dx(i,j), i=1, ncomp-1)
    print*, 'Froude', (froud(i), i=1, ncomp)
    print*, 'Bed', (z(i,j), i=1, ncomp)
    print*, 'newArea', (newArea(i,j), i=1, ncomp)
    print*, 'I2_corr', (ci2(i), i=1, ncomp)
    print*, 'Courant no', (courant(i), i=1, ncomp-1)
    print*, 'Maximum Courant no', maxCourant

    !
!10  format(f12.2 , <ncomp>f12.3)
10  format(f12.2 ,i6, 1200f13.3)

    call cpu_time( t2 )
    print '("Time = ",f10.3," seconds.")',t2 - t1
    !pause 202
end program mesh

function r_interpol(x,y,kk,xt)

    integer, intent(in) :: kk
    real, intent(in) :: xt, x(kk), y(kk)

    if (xt.lt.maxval(x) .and. xt.ge.minval(x)) then
        do k=1,kk-1
            if((x(k)-xt)*(x(k+1)-xt).le.0)then

                yt=(xt-x(k))/(x(k+1)-x(k))*(y(k+1)-y(k))+y(k)

                EXIT
            endif
        end do
    else
        print*, xt, ' is not within the limit'
        print*, 'maxval(x)= ', maxval(x), 'and minval(x)=', minval(x),'so',  xt, ' is not within the limit'
        print*, 'kk', kk
        print*, 'x', (x(i), i=1, kk)
        print*, 'y', (y(i), i=1, kk)
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

