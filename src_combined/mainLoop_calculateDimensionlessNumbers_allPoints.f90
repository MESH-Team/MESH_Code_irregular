subroutine call_dimensionless_numbers

    use constants_module
    use arrays_module
    use var_module
    use arrays_section_module

    implicit none

    integer(kind=4) :: i
    real(kind=4) :: wl_us, depth_us, q_us, v_us, pere_us, r_us, sk_us, ch_us
    real(kind=4) :: wl_ds, depth_ds, q_ds, v_ds, pere_ds, r_ds, sk_ds, ch_ds
    real(kind=4) :: ch_star_avg, channel_length, avg_celerity, avg_velocity, avg_depth
    real(kind=4) :: maxValue, dimlessWaveLength

    !avg_celerity = celerity(2)

    maxValue = 1e7

    do i=1, ncomp-1

        avg_celerity = (celerity2(i) + celerity2(i+1)) / 2.0    ! 'celerity2' is calculated celerity. 'celerity' is spatially averaged

        wl_us = newY(i)
        depth_us = newArea(i) / bo(i)
        q_us = newQ(i)
        v_us = abs( newQ(i) / newArea(i) )
        pere_us = pere(i)
        r_us = newArea(i) / pere(i)
        sk_us = sk(i)
        ch_us = sk(i) * r_us ** (1./6.)


        wl_ds = newY(i+1)
        depth_ds = newArea(i+1) / bo(i+1)
        q_ds = newQ(i+1)
        v_ds = abs( newQ(i+1) / newArea(i+1) )
        pere_ds = pere(i+1)
        r_ds = newArea(i+1) / pere(i+1)
        sk_ds = sk(i+1)
        ch_ds = sk(i+1) * r_ds ** (1./6.)


        ch_star_avg = ((ch_us + ch_ds) / 2.)  / sqrt( grav ) !! CORRECTED
        channel_length = dx(i)



        !wl_ds = newY(752)
        !depth_ds = newArea(752) / bo(752)
        !q_ds = newQ(752)
        !v_ds = abs( newQ(1) / newArea(752) )
        !pere_ds = pere(752)
        !r_ds = newArea(752) / pere(752)
        !sk_ds = sk(752)
        !ch_ds = sk(752) * r_ds ** (1./6.)


        !ch_star_avg = ((ch_us + ch_ds) / 2.)  / sqrt( grav ) !! CORRECTED
        !channel_length = sum(dx(1:752-1))

        dimlessWaveLength= 4000. !! new



        avg_velocity = (v_us + v_ds) / 2.
        avg_depth = (depth_us + depth_ds) / 2.

        dimensionless_Cr(i) = abs(avg_velocity / avg_celerity)
        if (dimensionless_Cr(i) .gt. maxValue) dimensionless_Cr(i) = maxValue

        dimensionless_Fo(i) = avg_velocity / sqrt(grav * avg_depth)
        if (dimensionless_Fo(i) .gt. maxValue) dimensionless_Fo(i) = maxValue

        !dimensionless_Fi(i) = 2*dimensionless_Cr(i) / (ch_star_avg ** 2.) * (channel_length / avg_depth) !! CORRECTED
        dimensionless_Fi(i) = 2*dimensionless_Cr(i) / (ch_star_avg ** 2.) * dimlessWaveLength !(channel_length / avg_depth) !! CORRECTED
        if (dimensionless_Fi(i) .gt. maxValue) dimensionless_Fi(i) = maxValue

        dimensionless_Fc(i) = dimensionless_Cr(i) * dimensionless_Fi(i)
        if (dimensionless_Fc(i) .gt. maxValue) dimensionless_Fc(i) = maxValue

        dimensionless_Di(i) = (dimensionless_Cr(i) / dimensionless_Fo(i)) ** 2. !! CORRECTED
        if (dimensionless_Di(i) .gt. maxValue) dimensionless_Di(i) = maxValue

        dimensionless_D(i)  = dimensionless_Di(i) / dimensionless_Fc(i)
        if (dimensionless_D(i) .gt. maxValue) dimensionless_D(i) = maxValue

        !print*, i, avg_velocity, avg_celerity, dimensionless_Cr(i)

    end do



    !print*, ch_star_avg, dimensionless_Cr, dimensionless_Fo, dimensionless_Fi
    !pause 5000

end subroutine
