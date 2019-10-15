!            FINITE DIFFERENCE METHOD
!
!  A program for one dimensional flow in open channel
!

program mesh

    use constants_module
    use arrays_module
    use arrays_section_module
    use var_module
    use matrix_module
    use sgate_module

    implicit none

    ! Local storage
    integer :: i, n, ntim, igate, igatmaxTableLengthe, pp, jj
    real(kind=4) :: areanew, cour, da, dq, dxini, yy, x, thetas, thesinv, tfin
    real(kind=4) :: t, skk, qq, qn, xt, r_interpol
    real(kind=4) :: arean, areac, hyrdn, hyrdc, perimn, perimc, qcrit, s0ds

    character(len=128) :: input_path, upstream_path , downstream_path
    character(len=128) :: bed_elevation_path, output_path
    character(len=128) :: channel_width_path
    character(len=128) :: path

    ! open file for input data

    !if (command_argument_count() > 0) then
    !    call get_command_argument(1, input_path)
    !    input_path = trim(input_path)
    !else
       ! input_path =arrays_module
    !end if
    ! new from Nazmul
    !character(len=128) :: all_path
    !character(len=128) :: output_area
    character(len=128) :: dx_path  !, output_q, output_wl

    open(unit=1,file="../3/input/input.txt",status='unknown')

    ! read data
    read(1,*) dtini
    read(1,*) dxini
    read(1,*) tfin
    read(1,*) ncomp
    read(1,*) ntim
    read(1,*) phi
    read(1,*) theta
    read(1,*) thetas
    read(1,*) thesinv
    read(1,*) alfa2
    read(1,*) alfa4
    read(1,*) f
    read(1,*) skk
    read(1,*) yy
    read(1,*) qq
    read(1,*) cfl
    read(1,*) ots
    read(1,*) yw
    read(1,*) bw
    read(1,*) w
    read(1,*) option
    read(1,*) yn
    read(1,*) qn
    read(1,*) igatmaxTableLengthe
    read(1,*) xSection_path
    read(1,*) bed_elevation_path
    read(1,*) upstream_path
    read(1,*) downstream_path
    read(1,*) channel_width_path
    read(1,*) dx_path
    read(1,*) output_path
    read(1,*) option_dsbc
    read(1,*) maxTableLength
    read(1,*) nel
    close (1)

    ! Allocate arrays
    call setup_arrays(ntim, ncomp)
    ! added Nazmul
    ! Allocate arrays
    call setup_arrays_section

    dt = dtini

    open(unit=90, file=trim(dx_path))
    do i=1,ncomp-1
        read(90, *) x, dx(i)
        !print*, x, dx(i)
        !pause 1
    end do
    close(90)
    !dx = dxini
    !	new addition Nazmul
    do i=1,ncomp
        call readXsection(i,(1.0/skk))
        ! Nazmul: This subroutine creates attribute table for each cross sections and saves in the hdd
        ! setting initial condition
        y(1,i) = z(i) + yy
        !print*, 'z=',z(i),'y(1,i)=',y(1,i)
        !pause 1
    end do

    ityp = 1

    ! setting initial condition
    q(1, :) = qq
    sk = skk
    x = 0.0

    ! Read hydrograph input Upstream
    open(unit=98, file=upstream_path)
    do n=1,ntim
      read(98,*) t, q(n, 1)
    end do
    close(98)
    ! Read hydrograph input downstream
    open(unit=99, file=downstream_path)
    do n=1,ntim
      read(99,*) t, y(n, ncomp)
    end do
    close(99)

    ! Read hydrograph input
    ! Nazmul: for irregular channel, bo = topwidth. It changes
    ! at every time step. Should be defined inside section(n)

    !open(unit=100, file=channel_width_path)
    !open(unit=100, file="./channel_width.txt")
    !do i=1, ncomp
    ! read(100,*) x, bo(i)
    !end do
    !close(100)


    ! Open files for output
    !output_path="../output/"
    !path = trim(output_path) // 'area.txt'
    path = trim(output_path) // 'output_wl.txt'
    open(unit=8, file=trim(path), status='unknown')
    path = trim(output_path) // 'q.txt'
    open(unit=9, file=trim(path), status='unknown')

    ! Output initial condition
    t = 0.0
    !write(8, *) t, ((y(1, i) - z(i)) * bo(i), i=1,ncomp)
    write(8, *) t, (y(1, i), i=1,ncomp)
    !print*, t, (y(1, i), i=1,ncomp)
    write(9, *) t, (q(1, i), i=1, ncomp)

    !
    ! Loop on time
    !
    do n=1, ntim-1

        ! Set upstream discharge
        dqp(1) = q(n+1,1) - q(n,1)

        call section(n)
        ! Nazmul: The subroutine calls the attribute tables and interpolate according to the available water level
        thes=thetas
        call matrixp(n)

        do i=2,ncomp
            cour=dt(i)/dx(i-1)
            rhs1=-cour*(f1(i)-f1(i-1)-d1(i)+d1(i-1))
            rhs2=-cour*(f2(i)-f2(i-1)-d2(i)+d2(i-1))+dt(i)*grav*(ci2(i)+aso(i))
            c11=g11inv(i)*b11(i-1)+g12inv(i)*b21(i-1)
            c12=g11inv(i)*b12(i-1)+g12inv(i)*b22(i-1)
            c21=g21inv(i)*b11(i-1)+g22inv(i)*b21(i-1)
            c22=g21inv(i)*b12(i-1)+g22inv(i)*b22(i-1)
            dap(i)=g11inv(i)*rhs1+g12inv(i)*rhs2-c11*dap(i-1)-c12*dqp(i-1)
            dqp(i)=g21inv(i)*rhs1+g22inv(i)*rhs2-c21*dap(i-1)-c22*dqp(i-1)
        end do

        ! Boundary conditions at downstream (right boundary)
        if (option_dsbc.eq.1) then
            dac(ncomp)=dap(ncomp)
            yn=(area(ncomp)+dap(ncomp))/bo(ncomp)
            arean=yn*bo(ncomp)
            perimn=2.0*yn+bo(ncomp)
            hyrdn=arean/perimn
            s0ds=-((z(ncomp)-z(ncomp-1))/dx(ncomp))
            qn=skk*arean*hyrdn**(2.0/3.0)*sqrt(s0ds)
            dqp(ncomp)=qn-q(n,ncomp)
            dqc(ncomp)=dqp(ncomp)

        elseif(option_dsbc.eq.2)then
            dac(ncomp)=dap(ncomp)
            yn=(area(ncomp)+dap(ncomp))/bo(ncomp)
            areac=yn*bo(ncomp)
            perimc=2.0*yn+bo(ncomp)
            hyrdc=areac/perimc
            s0ds=-((z(ncomp)-z(ncomp-1))/dx(ncomp))
            qn=skk*areac*hyrdc**(2.0/3.0)*sqrt(s0ds)
            qcrit=1.05*(((yn**3.0)*(bo(ncomp)**2.0)*grav)**(1.0/2.0))
            write(*,*)qcrit
            dqp(ncomp)=qcrit-q(n,ncomp)
            dqc(ncomp)=dqp(ncomp)

        else
            dac(ncomp)=0.0
            dap(ncomp)=0.0
            dqc(ncomp)=dqp(ncomp)

        endif


        ! Update via predictor
        areap = area + dap
        qp = q(n, :) + dqp

        call secpred()
        thes=thesinv
        call matrixc()

        do i=ncomp-1,1,-1
            cour=dt(i)/dx(i)
            rhs1=-cour*(f1(i+1)-f1(i)-d1(i+1)+d1(i))
            rhs2=-cour*(f2(i+1)-f2(i)-d2(i+1)+d2(i))+dt(i)*grav*(ci2(i)+aso(i))
            c11=g11inv(i)*b11(i+1)+g12inv(i)*b21(i+1)
            c12=g11inv(i)*b12(i+1)+g12inv(i)*b22(i+1)
            c21=g21inv(i)*b11(i+1)+g22inv(i)*b21(i+1)
            c22=g21inv(i)*b12(i+1)+g22inv(i)*b22(i+1)
            dac(i)=g11inv(i)*rhs1+g12inv(i)*rhs2-c11*dac(i+1)-c12*dqc(i+1)
            dqc(i)=g21inv(i)*rhs1+g22inv(i)*rhs2-c21*dac(i+1)-c22*dqc(i+1)
        end do

        ! Upstream boundary condition
        ! Prescribed discharge at the upstream
        ! Area correction is calculated
        dqc(1)=dqp(1)

        ! Final update
        do i=1,ncomp
            da=(dap(i)+dac(i))/2.0
            dq=(dqp(i)+dqc(i))/2.0
            areanew=da+area(i)
            if(areanew <= 0.0) areanew=0.0001

!       Nazmul: Now calculate y based on area calculated
!-------------------------------------
            write(file_num,'(i3.3)')i

            open(unit=29,file=trim(xSection_path)//file_num//'_tab')

            read(29,*)

            do pp=1,maxTableLength
                read(29,*,end=300) elevTable(pp),areaTable(pp)

            enddo
300         close(29)

            jj=pp-1

    !       interpolate the cross section attributes based on FINAL CALCULATED area
            xt=areanew

            y(n+1,i)=r_interpol(areaTable,elevTable,jj,xt)

!-------------------------------------
!        y(n+1,i)=areanew/bo(i)+z(i)

            q(n+1,i)=q(n,i)+dq
        end do

        t = t + dtini
        print "('- cycle',i6,'  terminated')", n
        write(8, *) t, (y(n+1, i), i=1,ncomp)
        write(9, *) t, (q(n+1, i), i=1,ncomp)
    end do
    ! End of time loop

    close(8)
    close(9)

end program mesh
