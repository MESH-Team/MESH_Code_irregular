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
    use subtools

    implicit none

    ! Local storage
    integer(kind=4) :: i, j, k, ppn, qqn, n, ntim, igate, pp, boundaryFileMaxEntry, saveFrequency
    integer(kind=4) :: linknb_ds, linknb_us
    real(kind=4) :: qnp1_ds, qnp1_us, qsum, y_ds

    real(kind=4) :: cour, da, dq, x, saveInterval, width
    real(kind=4) :: qn, xt, maxCourant, dtini_given, nodenb, linknb
    real(kind=4) :: frds, areasum, yk_ncomp, yav, areak_ncomp, areav, sumOldQ, currentQ, area_ds
    real(kind=4) :: arean, areac, hyrdn, hyrdc, perimn, perimc, qcrit, s0ds, timesDepth, latFlowValue

    real(kind=4) :: t, r_interpol_time, tfin, t1, t2, t0 !t0 start time
    !doubleprecision ::  t, r_interpol_time, tfin, t1, t2, t0 !t0 start time

    integer(kind=4) :: tableLength
    real(kind=4) :: area_0, width_0, errorY, hydR_0, q_sk_multi

    real(kind=4) :: r_interpo_nn

    character(len=128) :: output_path, other_input, ndep_path
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
    !open(unit=1,file="../Rectangular_Y_Channel/input/input_naturalChannel_exact.txt",status='unknown')
    !open(unit=1,file="../Rectangular_Y_Channel/input/test.txt",status='unknown')
    !open(unit=1,file="../Rectangular_Y_Channel/input/input_naturalChannel_tidal.txt",status='unknown')
    !open(unit=1,file="../NHD_Y_Channel/input/input_naturalChannel_exact.txt",status='unknown')
    !open(unit=1,file="../Multijunction_Network/input/input_multiJunction.txt",status='unknown')
    !open(unit=1,file="../Multijunction_Network/input/input_multiJunction_NHD.txt",status='unknown')
    !open(unit=1,file="../Multijunction_Network/input/input_multiJunction_NHD_mixedRouting.txt",status='unknown')
    !open(unit=1,file="../Multijunction_Network/input/input_multiJunction_NHD_mixedRouting_8channel.txt",status='unknown')
    !open(unit=1,file="../Multijunction_Network/input/input_multiJunction_NHD_mixedRouting_11channel.txt",status='unknown')
    !open(unit=1,file="../Multijunction_Network/input/input_multiJunction_NHD_mixedRouting_16channel.txt",status='unknown')
    !open(unit=1,file="../Multijunction_Network/Large_NHD_Geometry/input_NHD_mixedRouting_29channel.txt",status='unknown')
    !open(unit=1,file="../Synthetic_Network_Test/input",status='unknown')
    !open(unit=1,file="../Rectangular_Y_Channel/input/input_naturalChannel_exact.txt",status='unknown')
    !open(unit=1,file="../Multijunction_Network/Large_NHD_Geometry_temp/input_NHD_mixedRouting_29channel_temp.txt",status='unknown')
    !open(unit=1,file="../NHD_Y_Channel/input/input_naturalChannel_exact_20201012.txt",status='unknown')
    !open(unit=1,file="../NHD_Y_Channel/input/input_naturalChannel_test.txt",status='unknown')
    open(unit=1,file="../lateralFlow_test/input/input_crank_nicolson_test_lateralFlow.txt",status='unknown')
    !open(unit=1,file="../NHD_Y_Channel/input/input_naturalChannel_test_1chn.txt",status='unknown')



    print*, 'Reading input file'

    ! read data
    read(1,*) dtini_given     ! in seconds
    dtini = dtini_given; lastKnownDiffuDT = dtini_given         !; print*, dtini; pause 500
    read(1,*) dxini
    read(1,*) t0        ! in hours
    read(1,*) tfin      ! in hours
	ntim = floor( (tfin - t0) / dtini * 3600)

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
	allocate(ini_y(nlinks))
	allocate(ini_q(nlinks))
	read(1,*) (ini_y(i), i=1, nlinks)!; print*, yy
	read(1,*) (ini_q(i), i=1, nlinks)
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

    read(1,*) nupbnds                       ! No of u/s boundary data files
    allocate(upBoundTableEntry(nupbnds))
    allocate(upstream_path(nupbnds))
    do i=1,nupbnds
        read(1,*) upstream_path(i)
    end do

    read(1,*) ndnbnds                      ! No of d/s boundary data files
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

    ! Reading lateral flow data starts
    allocate(noLatFlow(nlinks))
	read(1,*) (noLatFlow(i), i=1, nlinks)

    allocate(latFlowLocations(maxval(noLatFlow),nlinks))   ! all the first nodes where a lateral flow starts
    allocate(latFlowType(maxval(noLatFlow),nlinks))        ! Lateral flow type: Type = 1 = time series; Type 2 = flow as a function of upstream flow
    allocate(latFlowXsecs(maxval(noLatFlow),nlinks))       ! no of x-secs at the downstream that the lateral flow is applied

    do j = 1,nlinks
        ncomp=nx1(j)
        if (noLatFlow(j) .gt. 0) then
            read(1,*) (latFlowLocations(i,j), i=1, noLatFlow(j))
            !do i=1,noLatFlow(j)
            !    if ((latFlowLocations(i,j)-1)*(latFlowLocations(i,j)-ncomp) .eq. 0) then
            !        print*, 'ERROR: Lateral flow cannot be applied at the boundaries'
            !        stop
            !    end if
            !end do
        else
            read(1,*)
        end if
    end do
    do j = 1,nlinks
        ncomp=nx1(j)
        if (noLatFlow(j) .gt. 0) then
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
        else
            read(1,*)
        end if
    end do

    do j = 1,nlinks
        ncomp = nx1(j)
        if (noLatFlow(j) .gt. 0) then
            read(1,*) (latFlowXsecs(i,j), i=1, noLatFlow(j))
        else
            read(1,*)
        end if
    end do
    ! Reading lateral flow data ends

    ! Reading Q-SK table data data starts
	allocate(noQSKtable(nlinks))                                ! how many tables are there in each river reach
	read(1,*) (noQSKtable(i), i=1, nlinks)

    allocate(eachQSKtableNodeRange(2,maxval(noQSKtable),nlinks))
    do j = 1,nlinks
        if (noQSKtable(j) .gt. 0) then
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
        else
            read(1,*)
            read(1,*)
        end if
    end do
    ! Reading Q-SK table data data ends

    read(1,*) ndep_path ! Reading the location of network file
!print*, ndep_path; pause
    close (1)       ! all input data read is finished

    ! Allocate arrays
    call setup_arrays(ntim, maxval(nx1), maxTableLength, boundaryFileMaxEntry, maxval(noLatFlow), maxval(noQSKtable), nlinks)
    call setup_arrays_section
    call setup_xsec_attribute_module(nel, maxval(nx1),nlinks)

    dt = dtini
    minDx = 1e6

    do j = 1,nlinks
        ncomp = nx1(j)
        open(unit=90, file=trim(dx_path(j)), status='unknown')
            do i=1,ncomp-1
                read(90, *) x, dx(i,j)
            end do
            if (minval(dx(1:ncomp-1,j)) .le. minDx) minDx=minval(dx(1:ncomp-1,j))
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
            oldY(i,j) = ini_y(j)  + z(i,j)
            oldQ(i,j) = ini_q(j)
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
    !oldQ = qq

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



    !!! This part is the code from DongHa for network and then modified by Nazmul
    open(unit=80, file=trim(ndep_path), status='unknown')
    ! this file contains the information for network connectivity.
    ! left column shows river reach number (j) and
    ! the right column shows  the number of links that are immediately upstream of link j
    ! skip one header line
    read(80, *)
    do j = 1,nlinks
        read(80, *) x, ndep(j)
    end do

    ! ++++ Y channel connectivity ++++++!
    !* the number of links that are immediately upstream of link j
    !ndep(1)=0; ndep(2)=0; ndep(3)=0; ndep(4)=0; ndep(5)=2; ndep(6)=3; ndep(7)=1; ndep(8)=1


    allocate(uslinks(maxval(ndep),nlinks))
    uslinks = NAnum
    ! Next part of the ndep file reads the link number of k_th link that is immediately upstream of link j
    ! The first column is the link number j, and the rest of the columns are the k_th upstream link of j
    ! skip one header line or blank line
    read(80, *)
    do j = 1,nlinks
        read(80,*) x, (uslinks(i,j), i=1, ndep(j))
    end do

   ! uslinks(1,1)=NAnum  !not available
   ! uslinks(1,2)=NAnum  !not available
   ! uslinks(1,3)=NAnum
   ! uslinks(1,4)=NAnum
   ! uslinks(1,5)=1; uslinks(2,5)=2
   ! uslinks(1,6)=3; uslinks(2,6)=4; uslinks(3,6)=5
   ! uslinks(1,7)=6
   ! uslinks(1,8)=7


    ! Next part of the ndep file reads the link number of k_th link that is immediately downstream of link j
    ! skip one header line or blank line
    dslink = 0
    read(80, *)
    do j = 1,nlinks
        read(80,*) x, dslink(j)
    end do
    !* link number that is immediately downstream of link j
    ! dslink(1)=5; dslink(2)=5; dslink(3)=6; dslink(4)=6; dslink(5)=6; dslink(6)=7; dslink(7)=8


    !++++ in case of one channel +++!
    !* the number of links that are immediately upstream of link j
!    ndep(1)=0

    !* link number of k_th link that is immediately upstream of link j
!    uslinks(1,1)=NAnum  !not available
    !* link number that is immediately downstream of link j
!    dslink(1)=NAnum



    ! Next part of the ndep file reads the boundary condition of j_th link
    ! the first coulmn is the link number j, the 2nd column indicates the u/s condition of j,
    ! and the 3rd column indicates the d/s condition of j
    !* when data at either upper or lower end of link j is available,
    !* instrdflag(j,1)=1 when water level data is known at the upper end of link j
    !* instrdflag(j,1)=2 when discharge data is known
    !* instrdflag(j,1)=3 when rating curve is known

    !* instrdflag(j,2)=1 when water level data is known at the lower end of link j
    !* instrdflag(j,2)=2 when discharge data is known
    !* instrdflag(j,2)=3 when rating curve is known
    !* Otherwise, instrdflag(j,1/2)=0


    ! skip one header line or blank line
    read(80, *)
    do j = 1,nlinks
        read(80, *) x, instrdflag(j,1), instrdflag(j,2)
    end do



   ! instrdflag(1,1)=2; instrdflag(1,2)=0    !*discharge data is known at the upper end of link 1
   ! instrdflag(2,1)=2; instrdflag(2,2)=0    !*discharge is known at the upper end of link 2
   ! instrdflag(3,1)=2; instrdflag(3,2)=0    !*discharge data is known at the upper end of link 3
   ! instrdflag(4,1)=2; instrdflag(4,2)=0    !*discharge data is known at the upper end of link 4
   ! instrdflag(5,1)=0; instrdflag(5,2)=0    !*discharge data is known at the upper end of link 5
   ! instrdflag(6,1)=0; instrdflag(6,2)=0    !*stage data is known at the lower end of link 6
   ! instrdflag(7,1)=0; instrdflag(7,2)=0    !*stage data is known at the lower end of link 7
   ! instrdflag(8,1)=0; instrdflag(8,2)=1    !*stage data is known at the lower end of link 8

    close(80)

    ! Nazmul: Need to create a connectivity table like the following:
    ! p=Junction sequence no, q= number of river reaches connected to that junction, (riverReachSequenceNoInThatJunction(i),i=1,q), (streamOrderOfEachRivers(i),i=1,q)

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
    ppn = 1; qqn = 1        ! ppn and qqn indicates the sequence no of the boundary data
    do j = 1, nlinks
        ncomp = nx1(j)

        if (instrdflag(j,1) .eq. 2) then
            ! interpolation of boundaries at the desired time step
            oldQ(1,j)=r_interpol_time(USBoundary(1, 1:upBoundTableEntry(ppn), ppn), &
                USBoundary(2, 1:upBoundTableEntry(ppn), ppn),upBoundTableEntry(ppn),t)
            ppn = ppn +1
        end if

        if (instrdflag(j,2) .eq. 1) then
            ! interpolation of boundaries at the desired time step
            oldY(ncomp,j)=r_interpol_time(DSBoundary(1, 1:downBoundTableEntry(qqn), qqn), &
                DSBoundary(2, 1:downBoundTableEntry(qqn), qqn),downBoundTableEntry(qqn),t)
            qqn = qqn +1
            ncompElevTable = xsec_tab(1,:,ncomp,j)
            ncompAreaTable = xsec_tab(2,:,ncomp,j)
            xt=oldY(ncomp,j)
            call r_interpol(ncompElevTable,ncompAreaTable,nel,xt,oldArea(ncomp,j))
			if (oldArea(ncomp,j) .eq. -9999) then
				print*, 'At j = ',j,', i = ',ncomp, 'time =',t, 'interpolation of oldArea(ncomp,j) was not possible'
				stop
			end if
        end if
    end do



 !   if (instrdflag(j,1) .eq. 2) then        !! I.e. for river 1 and 2
 !       ! interpolation of boundaries at the desired time step at upstream Q boundaries from given time series
 !       newQ(1,j)=r_interpol_time(USBoundary(1, 1:upBoundTableEntry(ppn), ppn), &
 !           USBoundary(2, 1:upBoundTableEntry(ppn), ppn),upBoundTableEntry(ppn),t+dtini/60.)
 !       ppn = ppn +1
 !   end if




    ! Applying initial condition of area
    do j=1, nlinks
        ncomp = nx1(j)
        do i=1,ncomp
            elevTable = xsec_tab(1,1:nel,i,j)
            areaTable = xsec_tab(2,1:nel,i,j)
            call r_interpol(elevTable,areaTable,nel,oldY(i,j),oldArea(i,j))
			if (oldArea(i,j) .eq. -9999) then
				print*, 'At j = ',j,', i = ',i, 'time =',t, 'interpolation of oldArea(i,j) was not possible'
				stop
			end if
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
    volRemain = -999
    do j=1, nlinks
        ncomp = nx1(j)
        do i=1,ncomp-1
            volRemain(i,j) = (oldArea(i,j)+oldArea(i+1,j))/2.0*dx(i,j)
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

    path = trim(output_path) // 'currentRoutingNormal.txt'
    open(unit=98, file=trim(path), status='unknown')

    path = trim(output_path) // 'routingNotChanged.txt'
    open(unit=99, file=trim(path), status='unknown')

    path = trim(output_path) // 'volRemain.txt'
    open(unit=991, file=trim(path), status='unknown')


	! Some essential initial parameters for Diffusive Wave
	theta = 1.0
	qpx = 0.

    width = 10. !!! average width of MS from BR to BC
    celerity = 1.0
    maxCelerity = 1.0
    diffusivity = 10.


    !!! setting initial values of dimensionless parameters
    !dimensionless_Cr, dimensionless_Fo, dimensionless_Fi, dimensionless_Fc, dimensionless_Di, dimensionless_D
    dimensionless_Fi = 10.1
    dimensionless_Fc = 10.1
    dimensionless_D  = 0.1
    currentROutingDiffusive = 1

    ! parameters for diffusive vs partial diffusive
    currentRoutingNormal = 1
    routingNotChanged = 0





        ! Output initial condition
    do j=1, nlinks
        ncomp = nx1(j)

        write(8, 10)  t*60.0, j, (oldY(i,j), i=1,maxval(nx1))
        write(9, 10)  t*60.0, j, (oldQ(i,j), i=1, maxval(nx1))
        write(51, 10) t*60.0, j, (oldArea(i,j), i=1, maxval(nx1))
        write(941, 10) t*60.0, j, (dimensionless_Cr(i,j), i=1, maxval(nx1)-1)
        write(942, 10) t*60.0, j, (dimensionless_Fo(i,j), i=1, maxval(nx1)-1)
        write(943, 10) t*60.0, j, (dimensionless_Fi(i,j), i=1, maxval(nx1)-1)
        write(944, 10) t*60.0, j, (dimensionless_Di(i,j), i=1, maxval(nx1)-1)
        write(945, 10) t*60.0, j, (dimensionless_Fc(i,j), i=1, maxval(nx1)-1)
        write(946, 10) t*60.0, j, (dimensionless_D(i,j), i=1, maxval(nx1)-1)
        write(95, 10) t*60.0, j, (celerity2(i), i=1, ncomp)
        write(96, 10) t*60.0, j, (diffusivity2(i), i=1, ncomp)
        write(97, *) t*60.0, j, currentROutingDiffusive(j)
        write(98, *) t*60.0, j, (currentRoutingNormal(i,j), i=1, maxval(nx1)-1)
        write(99, *) t*60.0, j, (routingNotChanged(i,j), i=1, maxval(nx1)-1)
        write(991, *) t*60.0, j, (volRemain(i,j), i=1, maxval(nx1)-1)

    end do


    frus2 = 9999.
    notSwitchRouting=0
    minNotSwitchRouting = 2         ! works between Dynamic and Diffusive switching
    minNotSwitchRouting2 = 0        ! works between full Diffusive and partial Diffusive switching
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

        !print*, 'running ', j

        ncomp = nx1(j)

        !+++-- Checking the dtini for possible diffusive wave model and applying it to the model.
        if (j .eq. 1) call calculateDT(t0, t,saveInterval, cfl, tfin, maxCelerity,dtini_given)
        !dtini = dtini_given

        !+++---------------------------------------------------------------------------------+
        !+ Run predictor step though all the links
        !+++---------------------------------------------------------------------------------+
        if (instrdflag(j,1) .eq. 2) then        !! I.e. for river 1 and 2
            ! interpolation of boundaries at the desired time step at upstream Q boundaries from given time series
            newQ(1,j)=r_interpol_time(USBoundary(1, 1:upBoundTableEntry(ppn), ppn), &
                USBoundary(2, 1:upBoundTableEntry(ppn), ppn),upBoundTableEntry(ppn),t+dtini/60.)
            ppn = ppn +1
        end if

        !print*, t+dtini, newQ(1)

        if (instrdflag(j,2) .eq. 1) then        !! I.e. for river 3
            ! interpolation of boundaries at the desired time step at downstream WL boundaries from given time series
            newY(ncomp,j)=r_interpol_time(DSBoundary(1, 1:downBoundTableEntry(qqn), qqn), &
                DSBoundary(2, 1:downBoundTableEntry(qqn), qqn),downBoundTableEntry(qqn),t+dtini/60.)

            ncompElevTable = xsec_tab(1,:,ncomp,j)
            ncompAreaTable = xsec_tab(2,:,ncomp,j)
            xt=newY(ncomp,j)
            call r_interpol(ncompElevTable,ncompAreaTable,nel,xt,newArea(ncomp,j))
			if (newArea(ncomp,j) .eq. -9999) then
				print*, 'At j = ',j,', i = ',i, 'time =',t, 'interpolation of newArea(ncomp,j) was not possible'
				stop
			end if

            qqn = qqn +1


        end if


    !!START+++++++ If the channel has boundary originated from a junction+++++++
        !+++----------------------------------------------------------------
        !+ Hand over water from upstream to downstream properly according
        !+ to the nature of link connections, i.e., serial or branching.
        !+ Refer to p.52,RM1_MESH
        !+++----------------------------------------------------------------
        if (ndep(j).gt.0) then  !     !* the number of links that are immediately upstream of link j. Example: j = 3
            !*total water areas at n+1 at the end nodes of upstream links that join link j
            areasum=0.0
            do k=1, ndep(j)
                linknb=uslinks(k,j); nodenb=nx1(linknb)
                areasum=areasum + oldArea(nodenb,linknb) + dap(nodenb,linknb)
                !print*, 'areasum=', areasum!; pause
            end do
            dqp(1,j)=0.0;
            yav=0.0
            sumOldQ = 0.0







            do k=1, ndep(j)
                linknb=uslinks(k,j); nodenb=nx1(linknb)
                !**dqp(1,j)
                dqp(1,j)=dqp(1,j)+dqp(nodenb,linknb)    !! Not right! If initial condition is not as sum of Q is conversed, it will be wrong
                sumOldQ=sumOldQ+oldQ(nodenb,linknb)
                !**dap(1,j)
                !*area at the end nod of link k at time n+1
                areak_ncomp = oldArea(nodenb,linknb) + dap(nodenb,linknb)
                !bok_ncomp=bo(nodenb,linknb) !* bottom width at the end of link k in the immed. upstream of link j.
                !sslp=traps(nodenb, linknb)  !* side slope of trapezoidal channel.

                elevTable = xsec_tab(1,:,nodenb,linknb)
                areaTable = xsec_tab(2,:,nodenb,linknb)
                !print*, 'k, areak_ncomp',k, areak_ncomp!; pause
                call r_interpol(areaTable,elevTable,nel,areak_ncomp,yk_ncomp)
                !* weighted average based on areas at the end nodes of upstream link ks
                yav = yav + (areak_ncomp/areasum)*yk_ncomp
                !print*, 'yav', yav
                !print*, 'yav',yav, 'yk_ncomp', yk_ncomp!; pause
            end do
            !print*, 'Check boundary', j, newQ(1,j), newY(ncomp,j); pause


            dqp(1,j)=dqp(1,j)+sumOldQ-oldQ(1,j) ! Change from DongHa
            newQ(1,j)=oldQ(1,j)+dqp(1,j)        ! If this data is needed in diffusive wave, it takes newQ
            !* Area estimated for time n+1
            elevTable = xsec_tab(1,:,1,j)
            areaTable = xsec_tab(2,:,1,j)
            call r_interpol(elevTable,areaTable,nel,yav,areav)
            dap(1,j) = areav - oldArea(1,j)
        else            ! There are no links at the upstream of the reach. Example: j = 1, 2
            ! Set upstream discharge
            dqp(1,j) = newQ(1,j) - oldQ(1,j)
            dap(1,j) = 0.0
        end if
      !!END+++++++ If the channel has boundary originated from a junction+++++++



          !print*, 'Check boundary', 'areav', areav

        ! new approach to add multiple lateral flow to the same node
        lateralFlow = 0
        do i=1,noLatFlow(j)
            if (latFlowType(i,j) .eq. 1) then
                latFlowValue = r_interpol_time(lateralFlowTable(1, 1:dataInEachLatFlow(i,j), i,j), &
                    lateralFlowTable(2, 1:dataInEachLatFlow(i,j), i,j),dataInEachLatFlow(i,j),t)
                !print*, t, latFlowValue
            elseif (latFlowType(i,j) .eq. 2) then
                call r_interpol(lateralFlowTable(1, 1:dataInEachLatFlow(i,j), i,j), &
                    lateralFlowTable(2, 1:dataInEachLatFlow(i,j), i,j),dataInEachLatFlow(i,j), &
                    oldQ(latFlowLocations(i,j)-1,j),latFlowValue)
            endif


            ! added condition for lateral flow at the upstream boundary node
            if (latFlowLocations(i,j) .eq. 1) then
                if (latFlowXsecs(i,j) .gt. 1) then
                    latFlowValue = latFlowValue / &
                    (dx(1,j)+sum(dx(latFlowLocations(i,j):latFlowLocations(i,j)+latFlowXsecs(i,j)-2,j)))
                    do k=1,latFlowXsecs(i,j)
                        lateralFlow(latFlowLocations(i,j)+k-1)=lateralFlow(latFlowLocations(i,j)+k-1) + latFlowValue
                    end do
                else
                    latFlowValue = latFlowValue / dx(1,j)
                    lateralFlow(1)=latFlowValue
                end if
                newQ(1,j) = newQ(1,j)+lateralFlow(1)*dx(1,j)
                dqp(1,j)  = newQ(1,j) - oldQ(1,j)

            else
!                if (j .eq. 1) print*, latFlowValue, sum(dx(latFlowLocations(i,j)-1:latFlowLocations(i,j)-1+latFlowXsecs(i,j)-1,j))
                latFlowValue = latFlowValue / &
                        !sum(dx(latFlowLocations(i):latFlowLocations(i)+latFlowXsecs(i)-1))          !!! check this line
                       sum(dx(latFlowLocations(i,j)-1:latFlowLocations(i,j)-1+latFlowXsecs(i,j)-1,j))  ! current line until 20201111
              !          (sum(dx(latFlowLocations(i,j)-1:latFlowLocations(i,j)-1+latFlowXsecs(i,j)-1,j)) &         ! change 20201111 ! did not work
              !          -dx(latFlowLocations(i,j)-1,j)/2.0+ dx(latFlowLocations(i,j)-1+latFlowXsecs(i,j),j)/2.0)  ! change 20201111 ! did not work
              !  print*, (sum(dx(latFlowLocations(i,j)-1:latFlowLocations(i,j)-1+latFlowXsecs(i,j)-1,j)) &
              !          -dx(latFlowLocations(i,j)-1,j)/2.0+ dx(latFlowLocations(i,j)-1+latFlowXsecs(i,j),j)/2.0)
                !pause

                !test
                latFlowValue = 0.01

                do k=1,latFlowXsecs(i,j)
                    lateralFlow(latFlowLocations(i,j)+k-1)=lateralFlow(latFlowLocations(i,j)+k-1) + latFlowValue
                end do
            end if
        end do






          !! Calculating Y normal and Y critical at the upstream!!
        S_0 = (-z(2,j)+z(1,j))/dx(1,j)
        if (S_0 .gt. 0.) then
            currentQ = newQ(1,j)
            ! Calculating the Q_sk_multiplier for thie currentQ
            call calc_q_sk_multi(1,j,currentQ,q_sk_multi)
            ! Calculating the normal depth and critical depth at the river reach upstream as an output
            call normal_crit_y(1,j,q_sk_multi, S_0, currentQ, y_norm_us, y_crit_us, area_norm, area_crit)
        end if
        !print*, j, 'y_norm_us',y_norm_us




        !print*, j, 'newQ(1,j)',newQ(1,j),lateralFlow!; pause










        ! checking the value of Fc and Fi in each river reach
        lowerLimitCount = 0; higherLimitCount = 0

        do i=1,ncomp-1
            if ((dimensionless_Fi(i,j) .ge. 5.) .or. (dimensionless_Fc(i,j)  .ge. 5.))then
                higherLimitCount(j) = higherLimitCount(j) + 1
            elseif ((dimensionless_Fi(i,j) .le. 3.) .or. (dimensionless_Fc(i,j)  .le. 3.))then
                lowerLimitCount(j) = lowerLimitCount(j) + 1
            end if
        end do


        ! new switching algorithm
        ! for now, auto switching of routing is disabled
        ! manual routing selection:
        ! higherLimitCount(j) = 0; lowerLimitCount(j) = ncomp         ! for dynamic
        ! higherLimitCount(j) = ncomp; lowerLimitCount(j) = ncomp         ! for diffusive


        !! Forcing all computation to diffusive routing
       ! higherLimitCount(1) = ncomp; lowerLimitCount(1) = ncomp
       ! higherLimitCount(2) = ncomp; lowerLimitCount(2) = ncomp
       ! higherLimitCount(3) = ncomp; lowerLimitCount(3) = ncomp
       ! higherLimitCount(4) = ncomp; lowerLimitCount(4) = ncomp
       ! higherLimitCount(5) = ncomp; lowerLimitCount(5) = ncomp
       ! higherLimitCount(6) = ncomp; lowerLimitCount(6) = ncomp
       ! higherLimitCount(7) = ncomp; lowerLimitCount(7) = ncomp
       ! higherLimitCount(8) = ncomp; lowerLimitCount(8) = ncomp
        higherLimitCount = ncomp; lowerLimitCount = ncomp



        !print*, 'check here',j, ncomp; pause
        ! running either dynamic or diffusive wave routing at each river reach
        if (higherLimitCount(j) .ge. ncomp/2.) then
            if ( (currentROutingDiffusive(j) .eq. 0) .and. (notSwitchRouting(j) .lt. minNotSwitchRouting)) then
                call mesh_dynamic_predictor(dtini_given, t0, t, tfin, saveInterval,j)
                currentROutingDiffusive(j) = 0
            else
                call mesh_diffusive_forward(dtini_given, t0, t, tfin, saveInterval,j)
                if (currentROutingDiffusive(j) .eq. 0) notSwitchRouting(j) = 0
                currentROutingDiffusive(j) = 1
            end if
        elseif (lowerLimitCount(j) .ge. ncomp/2.) then

            if ( (currentROutingDiffusive(j) .eq. 1) .and. (notSwitchRouting(j) .lt. minNotSwitchRouting)) then
                call mesh_diffusive_forward(dtini_given, t0, t, tfin, saveInterval,j)
                currentROutingDiffusive(j) = 1
            else
                call mesh_dynamic_predictor(dtini_given, t0, t, tfin, saveInterval,j)
                if (currentROutingDiffusive(j) .eq. 1) notSwitchRouting(j) = 0
                currentROutingDiffusive(j) = 0
            end if
        else
            if (currentROutingDiffusive(j) .eq. 1) then
                call mesh_diffusive_forward(dtini_given, t0, t, tfin, saveInterval,j)
                currentROutingDiffusive(j) = 1
            else
                call mesh_dynamic_predictor(dtini_given, t0, t, tfin, saveInterval,j)
                currentROutingDiffusive(j) = 0
            end if
        end if

        notSwitchRouting(j) = notSwitchRouting(j) + 1

    end do  ! end off j loop for predictor


!print*, 'qp',qp
!print*, 'ce',(celerity(i,2), i=1,ncomp)
!pause
     !!+++++++ corrector starts +++++++++++++

    do j =  nlinks,1,-1

        !print*, j


        ncomp = nx1(j)

!+++------------------------------------------------------------+
            !+ Handle downstream boundary condition for a link j that has a link
            !+ immediately downstream.
            !+ **Note. dac and dqc at the last node of the most downstream
            !+ link are computed in sub dsbc during predictor step.
            !+ Refer to p.53-1,RM1_MESH
            !+++------------------------------------------------------------+
            if (j<nlinks) then      ! i.e., the river reach is not the last reach. Example: j = 1 and j = 2
                linknb=dslink(j)    ! Which river reach is immediately downstream of reach j
                !*dac(ncomp,j)
                !print*, j,linknb; pause

                !area_ds=oldArea(1,linknb)+0.5*(dap(1,linknb)+dac(1,linknb))
                !xt = newArea(1,linknb)
                !elevTable = xsec_tab(1,:,1,linknb)
                !areaTable = xsec_tab(2,:,1,linknb)
                !call r_interpol(areaTable,elevTable,nel,xt,y_ds)
                newY(ncomp,j)= newY(1,linknb)

                xt = newY(ncomp,j)
                elevTable = xsec_tab(1,:,ncomp,j)
                areaTable = xsec_tab(2,:,ncomp,j)
                call r_interpol(elevTable,areaTable,nel,xt,newArea(ncomp,j))
                !dmy=area(1,linknb)+dac(1,linknb)

                dac(ncomp,j)=2*(newArea(ncomp,j)-oldArea(ncomp,j))-dap(ncomp,j)
                !* dqc(ncomp,j)
                dqc(ncomp,j)=dqc(1,linknb)*qp(ncomp,j)/qp(1,linknb) ! Changed from what DongHa originally proposed.

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

            else            !! The river reach is the last reach. i.e. j = 3
                !print*,'j=', j; pause 500
                dap(ncomp,j)=newArea(ncomp,j) - oldArea(ncomp,j)    !update from downstream time series
                dac(ncomp,j)=dap(ncomp,j)
                dqc(ncomp,j)=dqp(ncomp,j)
            end if



        if (currentROutingDiffusive(j) .eq. 0) then
            call mesh_dynamic_corrector(dtini_given, t0, t, tfin, saveInterval,j)
        elseif (currentROutingDiffusive(j) .eq. 1) then
            call mesh_diffusive_backward(dtini_given, t0, t, tfin, saveInterval,j)
        else
            print*, 'Something is wrong in reach ', j
            stop
        end if

        ! Checking maximum celerity at the end of one time loop
	maxCelerity = 0.
        if (j .eq. 1) then
            do i=1,nlinks
                maxCelerity = maxval(celerity(1:nx1(i),i))
            end do
        end if

        write(95, 10) t*60.+dtini,j, (celerity2(i), i=1, ncomp)
        write(96, 10) t*60.+dtini,j, (diffusivity2(i), i=1, ncomp)

        ! checking full vs partial diffusive; which routing is applied at which node
        write(98, *) t*60.+dtini, j, (currentRoutingNormal(i,j), i=1, maxval(nx1)-1)
        write(99, *) t*60.+dtini, j, (routingNotChanged(i,j), i=1, maxval(nx1)-1)

    end do  ! end of j loop

    do j = 1, nlinks
        ncomp = nx1(j)
        do i=1,ncomp
            froud(i)=abs(newQ(i,j))/sqrt(grav*newArea(i,j)**3.0/bo(i,j))
            if (i .lt. ncomp) then
                courant(i)=(newQ(i,j)+newQ(i+1,j))/(newArea(i,j)+newArea(i+1,j))*dtini/dx(i,j)
            endif
        enddo
        if (maxCourant .lt. maxval (courant(1:ncomp-1))) then
            maxCourant = maxval (courant(1:ncomp-1))
        endif
    enddo

    do j=1, nlinks
        ncomp = nx1(j)
        do i=1,ncomp-1
            volRemain(i,j) = (newArea(i,j)+newArea(i+1,j))/2.0*dx(i,j)
        end do
    end do



    t = t + dtini/60.
    !t = t0*60. + (n+1)*dtini/60.
    !print "('- cycle',i9,'  completed')", n
    !print "('- simulation time ',2f16.2,'  completed')", t


    !if(mod(n+1,24*saveFrequency) .eq. 0 .or. (n.eq.0))write(*,*)'Nstep =', n, 'Days = ', t/60./24., real(n+1)/real(ntim)*100., '% completed'
    !print*, 'dqp', (dqp(i), i=1, ncomp)
    !print*, 'dqc', (dqc(i), i=1, ncomp)
    !print*, 'qp', (qp(i), i=1, ncomp)
    !print*, 'oldQ', (oldQ(i), i=1, ncomp)
    !print*, 'newQ', (newQ(i), i=1, ncomp) celerity

    !if (mod(n+1,saveFrequency) .eq. 0 .or. n .eq. (ntim-1)) then    !!! Calculate Fc and Fi
    do j = 1, nlinks
        ncomp = nx1(j)
        call calc_dimensionless_numbers(j)
        if ( (mod( (t-t0*60.)*60.  ,saveInterval) .le. TOLERANCE) .or. ( t .eq. tfin *60. ) ) then
            write(941, 10) t*60.0,j, (dimensionless_Cr(i,j), i=1, maxval(nx1)-1)
            write(942, 10) t*60.0,j, (dimensionless_Fo(i,j), i=1, maxval(nx1)-1)
            write(943, 10) t*60.0,j, (dimensionless_Fi(i,j), i=1, maxval(nx1)-1)
            write(944, 10) t*60.0,j, (dimensionless_Di(i,j), i=1, maxval(nx1)-1)
            write(945, 10) t*60.0,j, (dimensionless_Fc(i,j), i=1, maxval(nx1)-1)
            write(946, 10) t*60.0,j, (dimensionless_D(i,j),  i=1, maxval(nx1)-1)
        end if
    enddo


    !print*, 'times',mod( (t-t0*60.)*60.  ,saveInterval), t, t0, dtini, (currentROutingDiffusive(j),j=1,nlinks)
    print "(A, f6.2, A)", 'Parcent completed:', (t-t0*60.)/(tfin*60.-t0*60.)*100, '%'
    if ( (mod( (t-t0*60.)*60.  ,saveInterval) .le. TOLERANCE) .or. ( t .eq. tfin *60. ) ) then

    do j = 1, nlinks
        ncomp = nx1(j)
        write(8, 10) t*60.,j, (newY(i,j), i=1,maxval(nx1))
        write(9, 10) t*60.,j, (newQ(i,j), i=1,maxval(nx1))
        write(51, 10) t*60.,j, (newArea(i,j), i=1, maxval(nx1))

        write(91, 10) t*60.,j, (bo(i,j), i=1, maxval(nx1))
        write(93, 10) t*60.,j, (pere(i,j), i=1, maxval(nx1))

        write(97, *) t*60.0,j, currentROutingDiffusive(j), notSwitchRouting(j)
        write(991, *) t*60.0, j, (volRemain(i,j), i=1, maxval(nx1)-1)
    end do

    end if

    ! update of Y, Q and Area vectors
    oldY   = newY
    newY=-999
    oldQ   = newQ
    newQ=-999
    oldArea= newArea
    newArea=-999
    bo=-999
    pere=-999




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
    close(97)
    close(98)
    close(99)
    close(991)

   print*, 'dx', (dx(i,1), i=1, ncomp-1)
   ! print*, 'Froude', (froud(i), i=1, ncomp)
   print*, 'Bed'
   do j = 1, nlinks
        ncomp=nx1(j)
        print*, (z(i,j), i=1, maxval(nx1))
   end do
   ! print*, 'newArea', (newArea(i,j), i=1, ncomp)
   ! print*, 'I2_corr', (ci2(i), i=1, ncomp)
   ! print*, 'Courant no', (courant(i), i=1, ncomp-1)
   ! print*, 'Maximum Courant no', maxCourant

    !
!10  format(f12.2 , <ncomp>f12.3)
10  format(f12.2 ,i6, 1200f13.3)

    call cpu_time( t2 )
    print '("Total CPU Time = ",f10.3," seconds.")',t2 - t1
    !print "(A, f10.3, A)", 'Total CPU Time = ', t2 - t1, ' seconds.'
    !pause 202
end program mesh



function r_interpol_time(x,y,jj,xt)

    integer(kind=4), intent(in) :: jj
    real(kind=4), intent(in) :: x(jj), y(jj)
    real(kind=4), intent(in) :: xt
    real(kind=4) :: yt
    !real(kind=8), intent(out) :: r_interpol_time


    if (xt.le.maxval(x) .and. xt.ge.minval(x)) then
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

