module arrays_module

    implicit none
    save

    real(kind=4), allocatable :: area(:), y(:, :), q(:, :), bo(:)
    real(kind=4), allocatable :: areap(:), qp(:), z(:), dqp(:)
    real(kind=4), allocatable :: av11(:), av12(:), av21(:), av22(:)
    real(kind=4), allocatable :: dqc(:), dap(:), dac(:), ci1(:), ci2(:)
    real(kind=4), allocatable :: aso(:), f1(:), f2(:), depth(:)
    real(kind=4), allocatable :: g11inv(:), g12inv(:), g21inv(:), g22inv(:)
    real(kind=4), allocatable :: b11(:), b12(:), b21(:), b22(:)
    real(kind=4), allocatable :: eps2(:), eps4(:), d1(:), d2(:), u(:), c(:)
    real(kind=4), allocatable :: sk(:), co(:), gso(:), dbdx(:)
    real(kind=4), allocatable :: dt(:), dx(:)

    real(kind=4), allocatable :: pere(:),dpda(:)

    integer, allocatable :: ityp(:)

contains

    ! Allocate storage for all of the arrays in this module based on the number
    ! of time steps and spatial points
    subroutine setup_arrays(num_time, num_points)

        implicit none

        ! Input
        integer, intent(in) :: num_time, num_points

        allocate(area(num_points))
        allocate(y(num_time, num_points))
        allocate(q(num_time, num_points))
        allocate(bo(num_points))

        allocate(pere(num_points))
        allocate(dpda(num_points))


        allocate(areap(num_points))
        allocate(qp(num_points))
        allocate(z(num_points))
        allocate(dqp(num_points))
        allocate(av11(num_points))
        allocate(av12(num_points))
        allocate(av21(num_points))
        allocate(av22(num_points))
        allocate(dqc(num_points))
        allocate(dap(num_points))
        allocate(dac(num_points))
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
        allocate(sk(num_points))
        allocate(co(num_points))
        allocate(gso(num_points))
        allocate(dbdx(num_points))
        allocate(dt(num_points))
        allocate(ityp(num_points))
        allocate(dx(num_points-1))


    end subroutine setup_arrays

end module arrays_module
