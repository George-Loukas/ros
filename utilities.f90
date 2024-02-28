module utilities
    use settings
    implicit none
    contains

    subroutine print_settings(config)
        use constants
        implicit none
        type(setting), intent(inout) :: config
        call print_separ_line()
        write (*,'(a)') "Printing ROS Configuration"
        write (*,'(a,I2)') "   Max Iter. : ", max_iterations
        call print_separ_line()
        write (*,'(a)') "Angles"
        write (*,'(a)') "   Degrees"
        write (*,'(a,f7.3)') "      Alpha : ", config%alpha
        write (*,'(a,f7.3)') "      Beta  : ", config%beta
        write (*,'(a,f7.3)') "      Gamma : ", config%gamma
        call deg2rad(config%alpha)
        call deg2rad(config%beta)
        call deg2rad(config%gamma)
        write (*,'(a)') "   Radians"
        write (*,'(a,f10.3)') "      Alpha : ", config%alpha
        write (*,'(a,f10.3)') "      Beta  : ", config%beta
        write (*,'(a,f10.3)') "      Gamma : ", config%gamma
        call print_separ_line()
        write (*,'(a,f10.3)') "Vinit    : ", config%vinit
        write (*,'(a,f10.8)') "Epsilon  : ", config%epsilon
        write (*,'(a,f10.3)') "DT       : ", config%dt
        call print_separ_line()
        write (*,'(a)') "Indices"
        write (*,'(a,I1)') "   Symmetry      : ", config%idx_symm
        write (*,'(a,I1)') "   Inches        : ", config%idx_inch
        write (*,'(a,I1)') "   Ground Effect : ", config%idx_ground_eff
        call print_separ_line()
        write (*,'(a)') "Mesh"
        write (*,'(a,I7)') "   # Points    : ", config%num_points
        write (*,'(a,I7)') "   # Panels    : ", config%num_panels
        write (*,'(a,I7)') "   # Triangles : ", config%num_triangles
        write (*,'(a,I7)') "   # Quads     : ", config%num_quads
        call print_separ_line()
    end subroutine print_settings

    subroutine deg2rad(angle)
        use constants
        implicit none
        real(8), intent(inout) :: angle
        angle  = angle * pi / 180.d0
    end subroutine deg2rad

    subroutine inch2meter(length)
        implicit none
        real(8), intent(inout) :: length
        length = length * 0.0254
    end subroutine inch2meter

    subroutine read_settings(filename, configuration)
        implicit none
        character(len=*), dimension(3), intent(in) :: filename
        type(setting), intent(out)   :: configuration
        
        open(unit = 10, file = filename(1), status = 'unknown')
        read(10,100) configuration%alpha,   &
        &            configuration%beta ,   &
        &            configuration%gamma,   &
        &            configuration%vinit,   &
        &            configuration%epsilon, &
        &            configuration%flight_height, &
        &            configuration%idx_symm, &
        &            configuration%idx_inch, &
        &            configuration%idx_ground_eff

        configuration%num_points = count_lines(filename(2))
        configuration%num_panels = count_lines(filename(3))
        close(10)
        100 FORMAT(F10.3,/,F10.3,/,F10.3,/,F10.3,/,F10.8,/,F10.3 &
        &          ,/,I10,/,I10,/,I10)
    end subroutine read_settings

    function count_lines(filename) result(nlines)
        implicit none
        character(len=*)    :: filename
        integer             :: nlines
        integer             :: io
        open(10,file=filename, iostat=io, status='old')
        if (io/=0) stop 'Cannot open file! '

        nlines = 0
        do
            read(10,*,iostat=io)
            if (io/=0) exit
            nlines = nlines + 1
        end do

        close(10)
    end function count_lines

    subroutine print_status(iteration, lift, drag, convergence)
        implicit none
        integer, intent(in) :: iteration
        real(8), intent(in) :: lift, drag, convergence

        call print_separ_line()
        write (*,'(a,i2)') "Iteration : ", iteration
        write (*,'(a,f10.3)') "Lift   : ", lift
        write (*,'(a,f10.3)') "Drag   : ", drag 
        write (*,'(a,f12.8)') "Error  : ", convergence
        call print_separ_line

    end subroutine print_status

    subroutine print_separ_line()
        write (*,'(a)') "--------------------------"
    end subroutine

    subroutine solve_axb(a, n, b)
        implicit none
        integer, intent(in) :: n
        real(8), dimension(n, n), intent(in)    :: a
        real(8), dimension(n)   , intent(inout) :: b
      
        integer :: info
        real(8), dimension(n) :: ipiv(n)
        ! Call DGESV to solve the system AX = B. X is saved in B
        call dgesv(n,1,a,n,ipiv,b,n,info)
        if (info .ne. 0) then
          print *, 'Error in dgesv: ', info
          ! Handle the error appropriately
        end if
      end subroutine solve_axb
end module utilities