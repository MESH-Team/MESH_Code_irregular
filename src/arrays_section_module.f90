module arrays_section_module

    implicit none
    save

    real(kind=4), allocatable :: elevTable(:),areaTable(:)
    real(kind=4), allocatable :: pereTable(:),rediTable(:)
    real(kind=4), allocatable :: convTable(:),topwTable(:)
    real(kind=4), allocatable :: nwi1Table(:),elevTable_1(:),topwTable_1(:),areaTable_1(:),dPdATable(:)
    real(kind=4), allocatable :: eleNewSection(:),eleNewSection_1(:),topWNewSection(:),topWNewSection_1(:)

    integer :: maxTableLength, nel

    character*3 :: file_num, file_num_1
    character(len=128) :: xSection_path, xSection_tab, xSection_tab_1

contains

    subroutine setup_arrays_section

        implicit none

        allocate(elevTable(maxTableLength))
        allocate(areaTable(maxTableLength))
        allocate(pereTable(maxTableLength))
        allocate(rediTable(maxTableLength))
        allocate(convTable(maxTableLength))
        allocate(topwTable(maxTableLength))
        allocate(nwi1Table(maxTableLength))
        allocate(dPdATable(maxTableLength))
        allocate(elevTable_1(maxTableLength))
        allocate(topwTable_1(maxTableLength))
        allocate(areaTable_1(maxTableLength))

        allocate(eleNewSection(nel))
        allocate(eleNewSection_1(nel))
        allocate(topWNewSection(nel))
        allocate(topWNewSection_1(nel))

    end subroutine setup_arrays_section
end module
