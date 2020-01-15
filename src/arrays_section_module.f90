module arrays_section_module

    implicit none
    save

    real(kind=4), allocatable :: elevTable(:),areaTable(:)
    real(kind=4), allocatable :: pereTable(:),rediTable(:)
    real(kind=4), allocatable :: convTable(:),topwTable(:)
    real(kind=4), allocatable :: nwi1Table(:),upstreamEleTable(:),upstreamTopwTable(:),dPdATable(:)
    real(kind=4), allocatable :: eleNewSection(:),eleNewSectionUpstream(:),topWNewSection(:),topWNewSection_upstream(:)
    real(kind=4), allocatable :: eleNewSectionDownstream(:),topWNewSection_downstream(:)
    real(kind=4), allocatable :: downstreamEleTable(:),downstreamTopwTable(:),downstreamAreaTable(:)

    real(kind=4), allocatable :: elevtable_1(:),topwTable_1(:),areaTable_1(:)
    real(kind=4), allocatable :: topWNewSection_1(:),eleNewSection_1(:)

! only in version 2 20191511
    real(kind=4), allocatable :: I2Tablep(:),I2Tablec(:)
    real(kind=4), allocatable :: upstreamI2Tablec(:), downstreamI2Tablep(:)
    real(kind=4), allocatable :: dbdxTablep(:), dbdxTablec(:)
    real(kind=4), allocatable :: currentSquareDepth(:), downstreamSquareDepth(:), upstreamSquareDepth(:)

    integer :: maxTableLength, nel

    character*4 :: file_num, file_num_1
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
        allocate(upstreamEleTable(maxTableLength))
        allocate(upstreamTopwTable(maxTableLength))

        allocate(eleNewSection(nel))
        allocate(eleNewSectionUpstream(nel))
        allocate(topWNewSection(nel))
        allocate(topWNewSection_upstream(nel))

        allocate(eleNewSectionDownstream(nel))
        allocate(topWNewSection_downstream(nel))
        allocate(downstreamEleTable(nel))
        allocate(downstreamTopwTable(nel))
        allocate(downstreamAreaTable(nel))

        allocate(elevtable_1(nel))
        allocate(topwTable_1(nel))
        allocate(areaTable_1(nel))
        allocate(topWNewSection_1(nel))
        allocate(eleNewSection_1(nel))

! only in version 2 20151115
        allocate(I2Tablep(nel))
        allocate(I2Tablec(nel))
        allocate(upstreamI2Tablec(nel))
        allocate(downstreamI2Tablep(nel))
        allocate(dbdxTablep(nel))
        allocate(dbdxTablec(nel))
        allocate(currentSquareDepth(nel))
        allocate(downstreamSquareDepth(nel))
        allocate(upstreamSquareDepth(nel))

    end subroutine setup_arrays_section
end module
