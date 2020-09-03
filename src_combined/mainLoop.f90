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
    integer(kind=4) :: i, j, k, n, ntim, igate, pp, boundaryFileMaxEntry, ppp, qqq, noLatFlow, noQSKtable, saveFrequency
    integer(kind=4) :: lowerLimitCount, higherLimitCount, minNotSwitchRouting

    real(kind=4) :: cour, da, dq, dxini, yy, x, thetas, thesinv, saveInterval, width
    real(kind=4) :: skk, qq, qn, xt, r_interpol, maxCourant, dtini_given
    real(kind=4) :: arean, areac, hyrdn, hyrdc, perimn, perimc, qcrit, s0ds, timesDepth, latFlowValue

    real(kind=4) :: t, r_interpol_time, tfin, t1, t2, t0 !t0 start time
    !doubleprecision ::  t, r_interpol_time, tfin, t1, t2, t0 !t0 start time

    integer(kind=4) :: tableLength, notSwitchRouting
    real(kind=4) :: area_0, width_0, errorY, hydR_0, q_sk_multi

    real(kind=4) :: r_interpo_nn
    integer(kind=4) :: currentROutingDiffusive

    character(len=128) :: upstream_path , downstream_path
    character(len=128) :: manning_strickler_path, output_path, other_input
    character(len=128) :: QSKtablePath, lateralFlow_path
    character(len=128) :: path

    !real(kind=4) :: pere(500)


    ! open file for input data
    character(len=128) :: dx_path

    call cpu_time( t1 )

    !open(unit=1,file="../lower_Mississippi/input/input_BR_2_BC_2009.txt",status='unknown')
    !open(unit=1,file="../lower_Mississippi/input/input_BR_2_BC_2009.txt",status='unknown')
    !open(unit=1,file="../lateralFlow_test/input/input_dynamic_lateralFlow.txt",status='unknown')
    !open(unit=1,file="../Mississippi_River_11_years_20200511/input/input_Mississippi_BR2SWP_dynamic_new.txt",status='unknown')
    !open(unit=1,file="../Vermelion_River/input/input_Vermelion_dynamic_20200526.txt",status='unknown')
    !open(unit=1,file="../../../MESH_code/4-A-1_US_2bound/input_naturalChannel.txt",status='unknown')
    !open(unit=1,file="../../../MESH_code/New_DSP_test3/input.txt",status='unknown')
    !open(unit=1,file="../Mississippi_River_11_years_20200511/input/input_Mississippi_BR2SWP_dynamic_new.txt",status='unknown')
    !open(unit=1,file="../Vermelion_River/input/input_Vermelion_dynamic_20200526.txt",status='unknown')
    !open(unit=1,file="../../../MESH_code/4-A-1_US_2bound/PureSuperCritical_Boundary_Switch_Test/input.txt",status='unknown')
    !open(unit=1,file="../../../MESH_code/4-A-1_US_2bound/PureSuperCritical_Boundary_Switch_Test2/input.txt",status='unknown')
    !open(unit=1,file="../NHDplus_at_US_lateralFlowOnly/input.txt",status='unknown')
    open(unit=1,file="../NHD_Y_Channel/input/input_naturalChannel_exact1_interpol.txt",status='unknown')

    print*, 'Reading input file'

    ! read data
    read(1,*) dtini_given     ! in seconds
    dtini = dtini_given; lastKnownDiffuDT = dtini_given       !; print*, dtini; pause 500
    read(1,*) dxini
    read(1,*) t0        ! in hours
    read(1,*) tfin      ! in hours
	ntim = floor( (tfin - t0) / dtini * 3600)
	print*, ntim
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
    read(1,*) saveInterval !; saveFrequency = saveInterval / dtini
    read(1,*) noLatFlow

    allocate(latFlowLocations(noLatFlow))   ! all the first nodes where a lateral flow starts
    allocate(latFlowType(noLatFlow))        ! Lateral flow type: Type = 1 = time series; Type 2 = flow as a function of upstream flow
    allocate(latFlowXsecs(noLatFlow))       ! no of x-secs at the downstream that the lateral flow is applied

    read(1,*) (latFlowLocations(i), i=1, noLatFlow)
    !do i=1,noLatFlow
        !if ((latFlowLocations(i)-1)*(latFlowLocations(i)-ncomp) .eq. 0) then
        !    print*, 'ERROR: Lateral flow cannot be applied at the boundaries'
        !    stop
        !end if
    !end do
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

    ! Allocate arrays
    call setup_arrays(ntim, ncomp, maxTableLength, boundaryFileMaxEntry, noLatFlow, noQSKtable)
    call setup_arrays_section
    call setup_xsec_attribute_module(nel, ncomp)

    dt = dtini

    open(unit=90, file=trim(dx_path), status='unknown')
    do i=1,ncomp-1
        read(90, *) x, dx(i)
    end do
    close(90)

    ! reading Strickler's coefficient at each section
    open(unit=85,file=trim(manning_strickler_path), status='unknown') !! //'Mannings_Stricklers_coeff.txt', status='unknown')
    do i=1,ncomp
        read(85, *) x, sk(i)
        call readXsection(i,(1.0/sk(i)),timesDepth)
        ! This subroutine creates attribute table for each cross sections and saves in the hdd
        ! setting initial condition
        !y(1,i) = yy ! + z(i)
        oldY(i) = yy  + z(i)
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
    !    read(91, *) oldQ(i), oldY(i), oldArea(i)
    !end do
    !close(91)

    !q(1, :) = qq
    oldQ = qq


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
        read(87,*,end=301) USBoundary(1, n), USBoundary(2, n)       !! time column is in hours
    end do
301 close(87)
    ppp = n-1   ! stores the number of entries in the boundary file

    ! Read hydrograph input Downstream
    open(unit=88, file=downstream_path)
    do n=1,boundaryFileMaxEntry
      read(88,*,end=302)  DSBoundary(1, n), DSBoundary(2, n)       !! time column is in hours
    end do
302 close(88)
    qqq = n-1   ! stores the number of entries in the boundary file

    t=t0*60.0     !! from now on, t is in minute

    ! applying boundary
    ! interpolation of boundaries at the initial time step
    oldQ(1)    =r_interpol_time(USBoundary(1, 1:ppp),USBoundary(2, 1:ppp),ppp,t)
    oldY(ncomp)=r_interpol_time(DSBoundary(1, 1:qqq),DSBoundary(2, 1:qqq),qqq,t)
    !print*, t, oldY(ncomp)

    ! DS Boundary treatment: from water level to area time series
    ncompElevTable = xsec_tab(1,:,ncomp)
    ncompAreaTable = xsec_tab(2,:,ncomp)

	!open(unit=81,file=trim(output_path)//'DS_area.txt', status='unknown')
    xt=oldY(ncomp)
    oldArea(ncomp)=r_interpol(ncompElevTable,ncompAreaTable,nel,xt)
    !write(81, *) t, oldArea(ncomp)

    ! Applying initial condition of area
    do i=1,ncomp-1
        elevTable = xsec_tab(1,:,i)
        areaTable = xsec_tab(2,:,i)
        oldArea(i)=r_interpol(elevTable,areaTable,nel,oldY(i))
    enddo


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

    path = trim(output_path) // 'froude.txt'
    open(unit=951, file=trim(path), status='unknown')

    path = trim(output_path) // 'diffusivity.txt'
    open(unit=96, file=trim(path), status='unknown')

    path = trim(output_path) // 'routing.txt'
    open(unit=97, file=trim(path), status='unknown')

    path = trim(output_path) // 'conveyence.txt'
    open(unit=98, file=trim(path), status='unknown')

    path = trim(output_path) // 'dAdX_dTdQ.txt'
    open(unit=99, file=trim(path), status='unknown')


	! Some essential initial parameters for Diffusive Wave
	theta = 1.0
	qpx = 0.

    width = 10. !!! average width of MS from BR to BC is 800; Average Vermilion width is 10
    celerity = 1.0
    diffusivity = 10.


    !!! setting initial values of dimensionless parameters
    !dimensionless_Cr, dimensionless_Fo, dimensionless_Fi, dimensionless_Fc, dimensionless_Di, dimensionless_D
    dimensionless_Fi = 10.1
    dimensionless_Fc = 10.1
    currentROutingDiffusive = 1

        ! Output initial condition

    write(8, 10)  t*60.0, (oldY(i), i=1,ncomp)
    write(9, 10)  t*60.0, (oldQ(i), i=1, ncomp)
    write(51, 10) t*60.0, (oldArea(i), i=1, ncomp)
    write(91, 10) t*60.0, (bo(i), i=1, ncomp)
    write(93, 10) t*60.0, (pere(i), i=1, ncomp)
    write(941, 10) t*60.0, (dimensionless_Cr(i), i=1, ncomp-1)
    write(942, 10) t*60.0, (dimensionless_Fo(i), i=1, ncomp-1)
    write(943, 10) t*60.0, (dimensionless_Fi(i), i=1, ncomp-1)
    write(944, 10) t*60.0, (dimensionless_Di(i), i=1, ncomp-1)
    write(945, 10) t*60.0, (dimensionless_Fc(i), i=1, ncomp-1)
    write(946, 10) t*60.0, (dimensionless_D(i), i=1, ncomp-1)
    write(95, 10) t*60.0, (celerity2(i), i=1, ncomp)
    write(951, 10) t*60., (froud(i), i=1, ncomp)
    write(96, 10) t*60.0, (diffusivity2(i), i=1, ncomp)
    write(97, *) t*60.0, currentROutingDiffusive, '0'
    write(98, 10) t*60.0, (co(i), i=1, ncomp)
    write(99, 10) t*60.0, (co(i)*0, i=1, ncomp)

    S_0 = (-z(2)+z(1))/dx(1)
    frus2 = 9999.
    notSwitchRouting=0
    minNotSwitchRouting = 0
    !
    ! Loop on time
    !
    !do n=0, ntim-1
    do while ( t .lt. tfin *60.)

        ! interpolation of boundaries at the desired time step
        newQ(1)     =r_interpol_time(USBoundary(1, 1:ppp),USBoundary(2, 1:ppp),ppp,t+dtini/60.)
    !print*, t+dtini, newQ(1)
        newY(ncomp) =r_interpol_time(DSBoundary(1, 1:qqq),DSBoundary(2, 1:qqq),qqq,t+dtini/60.)
    !print*, t+dtini, newY(ncomp)
        xt=newY(ncomp)
		newArea(ncomp)=r_interpol(ncompElevTable,ncompAreaTable,nel,xt)

		! new approach to add multiple lateral flow to the same node
        lateralFlow = 0
        do i=1,noLatFlow
            if (latFlowType(i) .eq. 1) then
                latFlowValue = r_interpol_time(lateralFlowTable(1, 1:dataInEachLatFlow(i), i), &
                    lateralFlowTable(2, 1:dataInEachLatFlow(i), i),dataInEachLatFlow(i),t)
                !print*, t, latFlowValue
            elseif (latFlowType(i) .eq. 2) then
                latFlowValue = r_interpol(lateralFlowTable(1, 1:dataInEachLatFlow(i), i), &
                    lateralFlowTable(2, 1:dataInEachLatFlow(i), i),dataInEachLatFlow(i),oldQ(latFlowLocations(i)-1))
            endif
            ! added condition for lateral flow at the upstream boundary
            if (latFlowLocations(i) .eq. 1) then
                latFlowValue = latFlowValue / &
                    (dx(1)+sum(dx(latFlowLocations(i):latFlowLocations(i)+latFlowXsecs(i)-1)))
                do k=1,latFlowXsecs(i)+1
                    lateralFlow(latFlowLocations(i)+k-1)=lateralFlow(latFlowLocations(i)+k-1) + latFlowValue
                end do
                newQ(1) = newQ(1)+lateralFlow(1)*dx(1)
            else
                latFlowValue = latFlowValue / &
                    sum(dx(latFlowLocations(i)-1:latFlowLocations(i)-1+latFlowXsecs(i)-1))

                do k=1,latFlowXsecs(i)
                    lateralFlow(latFlowLocations(i)+k-1)=lateralFlow(latFlowLocations(i)+k-1) + latFlowValue
                end do
            end if

        end do

        !print*, 'Lat Flow',lateralFlow


      !! Calculating Y normal at the upstream!!
      !temp disabled!
        if (S_0 .gt. 0.) then

            q_sk_multi = 1.0
              do pp = 1, size(Q_sk_tableEntry)
                if (  ( eachQSKtableNodeRange(1,pp) - 1) * ( eachQSKtableNodeRange(2,pp) - 1) .le. 0 ) then
                    tableLength = Q_sk_tableEntry(pp)
                    q_sk_multi = r_interpo_nn(Q_Sk_Table(1,1:tableLength,pp),Q_Sk_Table(2,1:tableLength,pp),tableLength,newQ(1))
                end if
              end do
              !slope = (z(1)-z(2))/dx(1)
              !call normal_crit_y(1, q_sk_multi, slope, newQ(1)), temp, temp2, newArea(i), temp3)


              elevTable = xsec_tab(1,:,1)
              areaTable = xsec_tab(2,:,1)
              rediTable = xsec_tab(4,:,1)
              topwTable = xsec_tab(6,:,1)
              area_0 = r_interpol(elevTable,areaTable,nel,oldY(1))
              width_0= r_interpol(elevTable,topwTable,nel,oldY(1))
              area_crit=area_0
              errorY = 100.
            do while (errorY .gt. 0.0001)

              hydR_0 = r_interpol(areaTable,rediTable,nel,area_0)
              area_norm = newQ(1)/sk(1)/q_sk_multi/ hydR_0 ** (2./3.) / sqrt(S_0)

              errorY = abs(area_norm - area_0) / area_0
              area_0 = area_norm
              area_crit= (newQ(1) * newQ(1) * width_0 / grav) ** (1./3.)
              width_0  = r_interpol(areaTable,topwTable,nel,area_crit)

            enddo

              !pause
          y_norm_us = r_interpol(areaTable,elevTable,nel,area_0)
          y_crit_us = r_interpol(areaTable,elevTable,nel,area_crit)
          print*, 'check point -1',newQ(1),1/sk(1),q_sk_multi, y_norm_us-z(1),  y_crit_us-z(1)
      endif



        lowerLimitCount = 0; higherLimitCount = 0

        do i=1,ncomp-1
            if ((dimensionless_Fi(i) .ge. 5.) .or. (dimensionless_Fc(i)  .ge. 5.))then
                higherLimitCount = higherLimitCount + 1
            elseif ((dimensionless_Fi(i) .le. 3.) .or. (dimensionless_Fc(i)  .le. 3.))then
                lowerLimitCount = lowerLimitCount + 1
            end if
        end do


        ! Testing of switch
        if (t .lt. 0) then
            if ( mod(floor((t - 60.)/600.),2) .lt. 1) then
                higherLimitCount = 0; lowerLimitCount = ncomp
            else
                higherLimitCount = ncomp; lowerLimitCount = ncomp
            end if
        else
            higherLimitCount = ncomp; lowerLimitCount = ncomp
        end if



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
        if (higherLimitCount .ge. ncomp/2.) then
            if ((currentROutingDiffusive .eq. 0) .and. (notSwitchRouting .lt. minNotSwitchRouting)) then
                call mesh_dynamic(dtini_given, ppp,qqq, t0, t, tfin, saveInterval)
                currentROutingDiffusive = 0
            else
                call mesh_diffusive(ppp,qqq, t0, t, tfin, saveInterval)
                if (currentROutingDiffusive .eq. 0) notSwitchRouting = 0
                currentROutingDiffusive = 1
            end if
        elseif (lowerLimitCount .ge. ncomp/2.) then
            if ((currentROutingDiffusive .eq. 1) .and. (notSwitchRouting .lt. minNotSwitchRouting)) then
                call mesh_diffusive(ppp,qqq, t0, t, tfin, saveInterval)
                currentROutingDiffusive = 1
            else
                call mesh_dynamic(dtini_given, ppp,qqq, t0, t, tfin, saveInterval)
                if (currentROutingDiffusive .eq. 1) notSwitchRouting = 0
                currentROutingDiffusive = 0
            end if
        else
            if (currentROutingDiffusive .eq. 1) then
                call mesh_diffusive(ppp,qqq, t0, t, tfin, saveInterval)
                currentROutingDiffusive = 1
            else
                call mesh_dynamic(dtini_given, ppp,qqq, t0, t, tfin, saveInterval)
                currentROutingDiffusive = 0
            end if
        end if

        !!! Calculate Fc and Fi
        call call_dimensionless_numbers


        do i=1,ncomp
            froud(i)=abs(newQ(i))/sqrt(grav*newArea(i)**3.0/bo(i))
            if (i .lt. ncomp) then
                courant(i)=(newQ(i)+newQ(i+1))/(newArea(i)+newArea(i+1))*dtini/dx(i)
            endif
        enddo

        if (maxCourant .lt. maxval (courant)) then
            maxCourant = maxval (courant)
        endif

        t = t + dtini/60.
        notSwitchRouting = notSwitchRouting + 1
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
        print*, 'times', t, t0, dtini, currentROutingDiffusive
        if ( (mod( (t-t0*60.)*60.  ,saveInterval) .eq. 0.0) .or. ( t .ge. tfin *60. ) ) then
            write(8, 10) t*60., (newY(i), i=1,ncomp)
            write(9, 10) t*60., (newQ(i), i=1,ncomp)
            write(51, 10) t*60., (newArea(i), i=1, ncomp)

            write(91, 10) t*60., (bo(i), i=1, ncomp)
            write(93, 10) t*60., (pere(i), i=1, ncomp)

            write(941, 10) t*60.0, (dimensionless_Cr(i), i=1, ncomp-1)
            write(942, 10) t*60.0, (dimensionless_Fo(i), i=1, ncomp-1)
            write(943, 10) t*60.0, (dimensionless_Fi(i), i=1, ncomp-1)
            write(944, 10) t*60.0, (dimensionless_Di(i), i=1, ncomp-1)
            write(945, 10) t*60.0, (dimensionless_Fc(i), i=1, ncomp-1)
            write(946, 10) t*60.0, (dimensionless_D(i),  i=1, ncomp-1)

            write(95, 10) t*60., (celerity2(i), i=1, ncomp)
            write(951, 10) t*60., (froud(i), i=1, ncomp)
            write(96, 10) t*60., (diffusivity2(i), i=1, ncomp)
            write(98, 10) t*60.0, (co(i), i=1, ncomp)
            write(99, 10) t*60.0, ((newArea(i)-oldArea(i))/dtini*dx(i-1)/(newQ(i)-oldQ(i)), i=2, ncomp)

        end if

        write(97, *) t*60.0, currentROutingDiffusive, notSwitchRouting
        !  if (n .eq. 8000) pause

        !write(81, *) t, newArea(ncomp)

        ! update of Y, Q and Area vectors
        oldY   = newY
        oldQ   = newQ
        oldArea= newArea

        !!! test
        !if (t .gt. 50000) dtini = 103
        !!! test end
    end do
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
    close(951)
    close(96)
    close(97)
    close(98)
    close(99)
    !close(81)

    print*, 'dx', (dx(i), i=1, ncomp-1)
    print*, 'Froude', (froud(i), i=1, ncomp)
    print*, 'Bed', (z(i), i=1, ncomp)
    print*, 'newArea', (newArea(i), i=1, ncomp)
    print*, 'I2_corr', (ci2(i), i=1, ncomp)
    print*, 'Courant no', (courant(i), i=1, ncomp-1)
    print*, 'Maximum Courant no', maxCourant

    !
!10  format(f12.2 , <ncomp>f12.3)
10  format(f14.3 , 1200f13.3)

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

