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



    !print*, 'a, b', a, b, time, initialTime, dtini, saveInterval ; pause 5000

    !if( (b .gt. a) .and. (time .gt. initialTime) ) then
    !    dtini = ( initialTime*60. + (a+1)*saveInterval/60. -  time ) * 60.
        !if (dtini .le. 1e-3) dtini = maxAllowCourantNo*minval(dx)/maxval(celerity)
    !end if

    if (b .gt. a) then
        dtini = (a+1) * ( saveInterval ) - (time-initialTime*60.)*60.
    end if

    if (dtini .gt. 100.) dtini = 100.

    if ( time+dtini/60. .gt. tfin*60. ) dtini =  (tfin*60.-time)*60.


    !print*, 'dt =', a, b,dtini!, time, saveInterval, maxAllowCourantNo

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

    if (b .gt. a) then
        dtini = (a+1) * ( saveInterval ) - (time-initialTime*60.)*60.
    end if

    if ( time+dtini/60. .gt. tfin*60. ) dtini =  (tfin*60.-time)*60.

end subroutine
