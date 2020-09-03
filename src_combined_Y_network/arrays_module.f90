module arrays_module

    implicit none
    save

    real, allocatable :: area(:), bo(:,:) !, y(:, :), q(:, :)
    !real, allocatable :: areap(:), qp(:), z(:), dqp(:)
    real, allocatable :: av11(:), av12(:), av21(:), av22(:)
    real, allocatable ::  ci1(:), ci2(:)!, dqc(:), dap(:), dac(:)
    real, allocatable :: aso(:), f1(:), f2(:), depth(:)
    real, allocatable :: g11inv(:), g12inv(:), g21inv(:), g22inv(:)
    real, allocatable :: b11(:), b12(:), b21(:), b22(:)
    real, allocatable :: eps2(:), eps4(:), d1(:), d2(:), u(:), c(:)
    real, allocatable :: co(:), gso(:), dbdx(:)!,sk(:)
    real, allocatable :: dx(:,:), froud(:), courant(:)
    real, allocatable :: dt(:)

	!**arrays for branching channel application
    integer, allocatable :: ndep(:), uslinks(:,:), dslink(:), instrdflag(:,:), nx1(:)
    !real, allocatable :: y(:, :, :), q(:, :, :), qlat(:,:,:), bo(:, :), traps(:,:) !,area(:,:), areafnal(:,:,:),
    real, allocatable :: areap(:, :), qp(:, :), z(:, :), sk(:, :)


    real, allocatable :: dqp(:,:), dap(:,:), dqc(:,:), dac(:,:)


    real, allocatable :: celerity(:,:), diffusivity(:,:), diffusivity2(:), celerity2(:)

    real, allocatable :: eei(:), ffi(:), exi(:), fxi(:), qpx(:,:), qcx(:)

    real, allocatable :: USBoundary(:,:,:), DSBoundary(:,:,:)
    integer, allocatable :: upBoundTableEntry(:), downBoundTableEntry(:)
! change for unsteady flow
    real, allocatable :: pere(:,:),dpda(:)

    real, allocatable :: oldQ(:,:), newQ(:,:), oldArea(:,:), newArea(:,:), oldY(:,:), newY(:,:)

    integer, allocatable :: ityp(:), latFlowLocations(:,:), dataInEachLatFlow(:,:), latFlowType(:,:), latFlowXsecs(:,:)
    real, allocatable :: lateralFlowTable(:,:,:,:), lateralFlow(:)


    real, allocatable :: dimensionless_Cr(:,:), dimensionless_Fo(:,:), dimensionless_Fi(:,:)
    real, allocatable :: dimensionless_Di(:,:), dimensionless_Fc(:,:), dimensionless_D(:,:)

    real, allocatable :: ini_y(:), ini_q(:)

    integer, allocatable :: Q_sk_tableEntry(:,:), noLatFlow(:), noQSKtable(:)
    real, allocatable :: eachQSKtableNodeRange(:,:,:), Q_sk_Table(:,:,:,:)

    real, allocatable :: lowerLimitCount(:), higherLimitCount(:)

    character(len=128), allocatable :: downstream_path(:), xSection_path(:), manning_strickler_path(:), upstream_path(:),dx_path(:)
    character(len=128), allocatable :: QSKtablePath(:), lateralFlow_path(:)

    integer, allocatable :: currentROutingDiffusive(:), notSwitchRouting(:)
    real :: minDx, maxCelerity

contains

    ! Allocate storage for all of the arrays in this module based on the number
    ! of time steps and spatial points
    subroutine setup_arrays(num_time, num_points, maxTableEntry1, maxTableEntry2, totalLatFlow, totalQSKtable, totalChannels)

        implicit none

        ! Input
        integer, intent(in) :: num_time, num_points, maxTableEntry1, maxTableEntry2, totalLatFlow, totalQSKtable, totalChannels

        allocate(area(num_points))

! change for unsteady flow

        allocate(bo(num_points,totalChannels))

        allocate(pere(num_points,totalChannels))
        allocate(dpda(num_points))


        allocate(areap(num_points,totalChannels))
        allocate(qp(num_points,totalChannels))
        allocate(z(num_points,totalChannels))
        allocate(av11(num_points))
        allocate(av12(num_points))
        allocate(av21(num_points))
        allocate(av22(num_points))

        allocate(dqp(num_points,totalChannels))
        allocate(dqc(num_points,totalChannels))
        allocate(dap(num_points,totalChannels))
        allocate(dac(num_points,totalChannels))

        allocate(ci1(num_points))
        allocate(ci2(num_points))
        allocate(aso(num_points))
        allocate(depth(num_points))
        allocate(f1(num_points))
        allocate(f2(num_points))
        allocate(g11inv(num_points))
        allocate(g12inv(num_points))
        allocate(g21inv(num_points))
        allocate(g22inv(num_points))
        allocate(b11(num_points))
        allocate(b12(num_points))
        allocate(b21(num_points))
        allocate(b22(num_points))
        allocate(eps2(num_points))
        allocate(eps4(num_points))
        allocate(d1(num_points))
        allocate(d2(num_points))
        allocate(u(num_points))
        allocate(c(num_points))
        allocate(sk(num_points,totalChannels))
        allocate(co(num_points))
        allocate(gso(num_points))
        allocate(dbdx(num_points))
        allocate(dt(num_points))
        allocate(ityp(num_points))
        allocate(dx(num_points-1,totalChannels))

        allocate(froud(num_points))


        allocate(Q_sk_Table(2, maxTableEntry1, totalQSKtable,totalChannels))
        allocate(Q_sk_tableEntry(totalQSKtable,totalChannels))


        allocate(USBoundary(2, maxTableEntry2,totalChannels))
        allocate(DSBoundary(2, maxTableEntry2,totalChannels))

        allocate(ndep(totalChannels))
        allocate(dslink(totalChannels))
        allocate(uslinks(num_points,totalChannels))

        allocate(instrdflag(totalChannels,2))

        allocate(courant(num_points-1))

        allocate(oldQ(num_points, totalChannels))
        allocate(newQ(num_points, totalChannels))
        allocate(oldArea(num_points, totalChannels))
        allocate(newArea(num_points, totalChannels))
        allocate(oldY(num_points, totalChannels))
        allocate(newY(num_points, totalChannels))

        allocate(lateralFlowTable(2, maxTableEntry2, totalLatFlow, totalChannels))
        allocate(dataInEachLatFlow(totalLatFlow, totalChannels))
        allocate(lateralFlow(num_points))



        allocate(celerity(num_points, totalChannels))
        allocate(diffusivity(num_points, totalChannels))
        allocate(celerity2(num_points))
        allocate(diffusivity2(num_points))


        allocate(eei(num_points))
        allocate(ffi(num_points))
        allocate(exi(num_points))
        allocate(fxi(num_points))
        allocate(qpx(num_points, totalChannels))
        allocate(qcx(num_points))

        allocate(dimensionless_Cr(num_points-1,totalChannels))
        allocate(dimensionless_Fo(num_points-1,totalChannels))
        allocate(dimensionless_Fi(num_points-1,totalChannels))
        allocate(dimensionless_Di(num_points-1,totalChannels))
        allocate(dimensionless_Fc(num_points-1,totalChannels))
        allocate(dimensionless_D(num_points-1,totalChannels))

        allocate(lowerLimitCount(totalChannels))
        allocate(higherLimitCount(totalChannels))

    end subroutine setup_arrays

end module arrays_module
