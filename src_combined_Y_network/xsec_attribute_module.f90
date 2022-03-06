module xsec_attribute_module

    implicit none
    save

    real, allocatable :: xsec_tab(:, :, :, :)

    contains

    ! Allocate storage for all of the arrays in this module based on the number
    ! of time steps and spatial points
    subroutine setup_xsec_attribute_module(elements, num_points, num_reaches)
        implicit none

        ! Input
        integer, intent(in) :: elements, num_points, num_reaches

        allocate(xsec_tab(12, elements, num_points, num_reaches))


        ! Headers:
        ! Elev(m)    Area(m2)     Peri(m)      Radi(m)   Conv(m3/s)    topWidth(m)    newI1(m3)    dPdA(1/m)     sk    dkdA    dbdxp(-)    dbdxc(-)      I2p     I2c



    end subroutine
end module xsec_attribute_module
