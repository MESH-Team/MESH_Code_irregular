subroutine calculateDT(initialTime, time, saveInterval, maxAllowCourantNo, tfin, max_C, given_dt)

    use var_module
    use arrays_module

    implicit none

    real, intent(in) :: initialTime, time, saveInterval, tfin, given_dt
    real, intent(in) :: maxAllowCourantNo, max_C
    integer          :: a, b

    !! initialTime is in hours
    !! tfin is in hours
    !! time is in minutes
    !! dtini is in seconds
    !! saveInterval is in seconds

    !print*, time, saveInterval, maxAllowCourantNo

    dtini = maxAllowCourantNo*minDx/max_C
    !dtini = min(dtini, given_dt)
    !dtini = given_dt

    a = floor( (time-initialTime*60.) /( saveInterval/60. ))            ! units:: time : min;  ! initialTime : hour ! saveInterval : sec
    b = floor( (time-initialTime*60.) +dtini/60.)/( saveInterval/60. )



    !print*, 'a, b', a, b, time, initialTime, dtini, saveInterval ; pause 5000

    if (b .gt. a) then
        dtini = (a+1) * ( saveInterval ) - (time-initialTime*60.)*60.
    end if

    if ( time+dtini/60. .gt. tfin*60. ) dtini =  (tfin*60.-time)*60.


    !! making dtini as a multiplier of 60
    !dtini = floor(dtini/60.) * 60

end subroutine


subroutine correctDT(initialTime, time, saveInterval, tfin)

    use var_module
    use arrays_module

    implicit none

    real(kind=4), intent(in) :: initialTime, time, saveInterval, tfin
    integer(kind=4)  :: a, b

    !! initialTime is in hours
    !! tfin is in hours
    !! time is in minutes
    !! dtini is in seconds
    !! saveInterval is in seconds

    !print*, time, saveInterval, maxAllowCourantNo

    a = floor( (time-initialTime*60.) /( saveInterval/60. ))
    b = floor(( (time-initialTime*60.) +dtini/60.)/( saveInterval/60. ))

    if( b .gt. a ) then
        dtini = ( initialTime*60. + (a+1)*saveInterval/60. -  time ) * 60.
    end if
    if ( a*saveInterval / 60. +initialTime*60.+dtini / 60. .gt. tfin*60. ) then
        dtini = ( tfin*60. - time )*60.
    end if

end subroutine
