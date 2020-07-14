subroutine calculateDT(initialTime, time, saveInterval, maxAllowCourantNo, tfin)

    use var_module
    use arrays_module

    implicit none

    real(kind=4), intent(in) :: initialTime, time, saveInterval, tfin
    real(kind=4), intent(in) :: maxAllowCourantNo
    integer(kind=4)  :: a, b

    !! initialTime is in hours
    !! tfin is in hours
    !! time is in minutes
    !! dtini is in seconds
    !! saveInterval is in seconds

    !print*, time, saveInterval, maxAllowCourantNo

    dtini = maxAllowCourantNo*minval(dx)/maxval(celerity)

    a = floor( (time-initialTime*60.) /( saveInterval/60. ))
    b = floor(( (time-initialTime*60.) +dtini/60.)/( saveInterval/60. ))

    if( b .gt. a ) then
        dtini = ( initialTime*60. + (a+1)*saveInterval/60. -  time ) * 60.
    end if
    if ( a*saveInterval / 60. +initialTime*60.+dtini / 60. .gt. tfin*60. ) then
        dtini = ( tfin*60. - time )*60.
    end if
    !print*, 'dt =', dtini, time, saveInterval, maxAllowCourantNo

end subroutine
