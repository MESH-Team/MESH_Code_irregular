module var_module

    implicit none
    save

    integer :: ncomp, ots, option_dsbc
    real(kind=4) :: cfl, qus, f, yds, rhs1, rhs2, c11, c12, c21, c22, us, thes
    real(kind=4) :: vv
    real(kind=4) :: dtini, lastKnownDiffuDT

    real(kind=4) :: frus2, S_0, y_norm_us, y_crit_us, area_norm, area_crit
end module var_module
