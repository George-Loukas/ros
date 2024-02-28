program main
    use constants
    use settings
    use utilities
    use geometry
    implicit none
    type(setting) :: config
    type(panel_) , allocatable :: panel(:)
    integer :: iteration = 0, npw = 0
    real    :: Tstart, Tend
    real(8) :: lift = 0.d0, prev_lift = 0.d0, drag = 0.d0, error = 1.d0
    integer :: idx_conv = 0
    real(8), allocatable, dimension(:,:) :: A
    real(8), allocatable, dimension(:)   :: B,C
    real(8), allocatable, dimension(:)   :: X1X,Y1Y,Z1Z

    real(8), allocatable, dimension(:,:) :: x_wake,y_wake,z_wake
    integer, allocatable, dimension(:,:) :: cw, icw
    integer, allocatable, dimension(:)   :: markw, ipropw
    ! Starting Timer
    call cpu_time(Tstart)

    call read_settings(filenames, config)
    ! config%num_panels = 172.3*1000 max memory 28/02/24 until ipropw allocation
    write (*,'(a)') 'Attempting to Allocate Memory'
    ! Allocating the geometry
    allocate(panel(config%num_panels))
    ! Allocating A B matrices
    allocate(A(config%num_panels,config%num_panels))
    allocate(B(config%num_panels))
    allocate(C(config%num_panels))
    ! Allocating the Wake
    allocate(x_wake(max_iterations,config%num_panels))   ! Maybe it doesnt need to be that big. Only the selected points from trailing edge
    allocate(y_wake(max_iterations,config%num_panels))   ! Maybe it doesnt need to be that big. Only the selected points from trailing edge
    allocate(z_wake(max_iterations,config%num_panels))   ! Maybe it doesnt need to be that big. Only the selected points from trailing edge
    allocate(cw(max_iterations,config%num_panels))       ! Maybe it doesnt need to be that big. Only the selected points from trailing edge
    allocate(icw(config%num_panels,2))                   ! Maybe it doesnt need to be that big. Only the selected points from trailing edge
    allocate(markw(config%num_panels))                   ! Maybe it doesnt need to be that big. Only the selected points from trailing edge
    allocate(ipropw(config%num_panels))                  ! Maybe it doesnt need to be that big. Only the selected points from trailing edge
    write (*,'(a)') 'Allocation Success'


    !call geometry_init()
    call print_settings(config)

    ! Zero Iteration - Preparation for main Loop
    !call vorcalc()
    
    ! First iteration
    iteration = 1
    call print_update_info()
    !call wake
    !CALL NEIBORG(NPAN,NGW)
    !CALL WAKREL(ITER,NPAN,NGW,NPW,DT,NSYM,NGRND)
    !CALL WAKINT(NPW,ITER)
    !call vorcalc()

    prev_lift = lift
    !call airload()
    call print_update_info()
    call check_convergence()
    iteration = iteration + 1

    do while (iteration<=max_iterations.and. idx_conv /= 1)
        !call wake
        !CALL WAKINT(NPW,ITER)
        !CALL WAKCOR(NPAN,NGW,ITER)
        !cALL WAKREL(ITER,NPAN,NGW,NPW,DT,NSYM,NGRND)
        !CALL WAKINT(NPW,ITER)
        !CALL WAKCOR(NPAN,NGW,ITER)
        !call vorcalc()
        prev_lift = lift
        !call airload()
        call print_update_info()
        call check_convergence()
        iteration = iteration + 1
    end do

    ! Stop Timer/Print Results
    call cpu_time(Tend)
    write (*,'(a,f10.2,a)') 'Execution time : ', Tend - Tstart, 's'

    contains

    subroutine print_update_info()
        implicit none
        logical :: info
        info = isnan(lift).or.lift==0
        info = info.and.(isnan(drag).or.drag==0)
        info = info.and.(isnan(error).or.error>1.d0)
        call print_separ_line()
        write(*,'(a,i2,a,i2)') 'Iteration : ', iteration,' of ',max_iterations
        if (info) then
            write(*,'(a,i1)') 'Lift  : ', 0
            write(*,'(a,i1)') 'Drag  : ', 0
            write(*,'(a,i1)') 'Error : ', 1
        else
            write(*,'(a,f10.3)') 'Lift  : ', lift
            write(*,'(a,f10.3)') 'Drag  : ', drag
            write(*,'(a,f10.5)') 'Error : ', error
        end if
        
        call print_separ_line()
    end subroutine print_update_info

    subroutine check_convergence()
        implicit none
        error = abs(prev_lift-lift)/lift
        if (error<=config%epsilon) idx_conv = 1
    end subroutine check_convergence

    subroutine geometry_init()
        implicit none
        real(8), allocatable, dimension(:,:) :: xyz
        integer, allocatable, dimension(:)   :: mark
        character(len = 30) :: format
        integer :: i,j,k,temp1,temp2,temp3,temp4, triangle_sum

        open(unit = 20, file = filenames(2), status = 'unknown')
        open(unit = 30, file = filenames(3), status = 'unknown')
        allocate(xyz(config%num_points,3), mark(config%num_points))
    
        format = '(I10,F15.5,F15.5,F15.5,I10)'
        do i = 1, config%num_points
            read (20,format) temp1,xyz(i,1),xyz(i,2),xyz(i,3),mark(i)
            if (config%idx_inch==1) then
                call inch2meter(xyz(i,1))
                call inch2meter(xyz(i,2))
                call inch2meter(xyz(i,3))
            end if
            config%dt = 0.25*abs(minval(xyz(:,1)) - maxval(xyz(:,1)))/config%vinit
        end do

        format = '(6I10)'
        triangle_sum = 0    ! counts triangles
        do i = 1, config%num_panels
            
            read (30,format) j, temp1,temp2,temp3,temp4, panel(i)%section
            
            k = abs(temp1)
            panel(i)%point(1,1) = xyz(k,1); panel(i)%point(1,2) = xyz(k,2); panel(i)%point(1,3) = xyz(k,3);
            if (mark(k) == 1) panel(i)%mark(1) = 1
            k = abs(temp2)
            panel(i)%point(2,1) = xyz(k,1); panel(i)%point(2,2) = xyz(k,2); panel(i)%point(2,3) = xyz(k,3);
            if (mark(k) == 1) panel(i)%mark(2) = 1
            k = abs(temp3)
            panel(i)%point(3,1) = xyz(k,1); panel(i)%point(3,2) = xyz(k,2); panel(i)%point(3,3) = xyz(k,3);
            if (mark(k) == 1) panel(i)%mark(3) = 1
            k = abs(temp4)
            panel(i)%point(4,1) = xyz(k,1); panel(i)%point(4,2) = xyz(k,2); panel(i)%point(4,3) = xyz(k,3);
            if (mark(k) == 1) panel(i)%mark(4) = 1
            if (temp3 == temp4) panel(i)%is_triangle = .true.
            if (panel(i)%is_triangle) triangle_sum = triangle_sum + 1

            call panel(i)%move_panel(config%beta,config%alpha,config%gamma,config%flight_height,config%idx_ground_eff)
            call panel(i)%get_centroid()
            call panel(i)%get_vectors()
        end do
        ! Counts all the triangles of the mesh and updates the config
        config%num_triangles = triangle_sum
        config%num_quads     = config%num_panels - triangle_sum
        close(20)
        close(30)
        deallocate(xyz,mark)
    end subroutine geometry_init

    subroutine vorcalc()
        use utilities
        implicit none        
        integer :: i = 0, j = 0, KPV = 0, NG = 0, NA = 0
        real(8) :: vxv = 0.d0, vyv = 0.d0, vzv = 0.d0
        real(8) :: vvx = 0.d0, vvy = 0.d0, vvz = 0.d0
        real(8) :: invz = 1.d0, invy= 1.d0
        real(8) :: VAIP = 0.d0

        
        if (iteration > 0) then ! not the first iteration
            B = 0.d0 ! resets B
            do i = 1, config%num_panels

                vxv=0.d0; vyv = 0.d0; vzv = 0.d0
                KPV = 1
                call velwak(panel(i)%centroid(1), &
                &           panel(i)%centroid(2), &
                &           panel(i)%centroid(3), &
                &           vxv,vyv,vzv)

                B(i) = - (config%vinit*KPV+vxv) * panel(i)%normal(1) &
                    &      -               vyv  * panel(i)%normal(2) &
                    &      -               vzv  * panel(i)%normal(3) &
                    &      + VAIP

            end do
        else
            A = 0.d0 ; B = 0.d0
            do i = 1, config%num_panels
                NG = 0; invz = 1.d0 ! if NG =0 we do not invert the Z (for ground effect)
                do while (.true.) 
                    do j = 1, config%num_panels
                        NA=0; invy = 1.d0 ! if NA = 0 we do not invert the Y (for symmetry)
                        do while (.true.)
                            call vortex(panel(i)%centroid(1), &
                            &           panel(i)%centroid(2), &
                            &           panel(i)%centroid(3), &
                            &           panel(j)%point(:,1),  &
                            &      invy*panel(j)%point(:,2),  &
                            &      invz*panel(j)%point(:,3),  &
                            & vvx,vvy,vvz)

                            A(i,j) = A(i,j) + (vvx*panel(i)%normal(1)  &
                            &               +  vvy*panel(i)%normal(2)  &
                            &               +  vvz*panel(i)%normal(3)) &
                            &               *((-1.d0)**NA)*((-1.d0)**NG)

                            if ((config%idx_symm == 0).or.(NA == 1)) exit ! go to the next j panel
                            NA = 1; invy = - invy ! When NA is 1 we invert the Y (symmetry)
                        end do
                    end do
                    if ((config%idx_ground_eff == 0).or.(NG == 1)) exit
                    NG = 1; invz = - invz ! When NG = 1 we invert the Z (for ground effect)
                end do
                B(i) = - (config%vinit) * panel(i)%normal(1) + VAIP
            end do
        end if
        C = B
        call solve_axb(A,config%num_panels,C)
    end subroutine vorcalc

    ! Subroutine to calculate velocity components (VX, VY, VZ) due to a vortex ring
    subroutine vortex(gx, gy, gz, xx, yy, zz, vx, vy, vz)
        use constants
        implicit none
        
      
        ! Input arguments
        real(8), intent(in) :: GX, GY, GZ ! Coordinates of the vortex center
        real(8), intent(in) :: XX(4), YY(4), ZZ(4) ! Coordinates of four points defining the vortex ring
        real(8), intent(inout) :: vx, vy, vz ! Output velocity components
      
        ! Local variables
        real(8) :: ABX, ABY, ABZ, AB  ! Components and magnitude of vector A
        real(8) :: APX, APY, APZ, AP  ! Components and magnitude of vector AP
        real(8) :: BPX, BPY, BPZ, BP  ! Components and magnitude of vector BP
        real(8) :: VPX, VPY, VPZ, VP  ! Components and magnitude of velocity
        real(8) :: ABPX, ABPY, ABPZ, ABP ! Components and magnitude of vector ABP
        real(8) :: V1, V2, V3  ! Intermediate variables for velocity calculation
        real(8) :: H           ! Norm of area element
        real(8) :: COSTH1, COSTH2 ! Cosine of angles

        integer :: i,j
      
        ! Initialize velocity components to zero
        VX = 0.d0
        VY = 0.d0
        VZ = 0.d0
      
        ! Loop over four points defining the vortex ring
        do i = 1, 4
            j = i + 1
            if (i == 4) j = 1  ! Handle periodic boundary condition

            ! Calculate vector A (difference between points j and i)
            ABX = XX(j) - XX(i)
            ABY = YY(j) - YY(i)
            ABZ = ZZ(j) - ZZ(i)
            AB = sqrt(ABX**2 + ABY**2 + ABZ**2)

            ! Calculate vector AP (difference between vortex center and point i)
            APX = GX - XX(i)
            APY = GY - YY(i)
            APZ = GZ - ZZ(i)
            AP = sqrt(APX**2 + APY**2 + APZ**2)

            ! Calculate COSTH1 (cosine of angle between A and AP)
            COSTH1 = (ABX * APX + ABY * APY + ABZ * APZ) / (AB * AP)

            ! Calculate vector BP (difference between vortex center and point j)
            BPX = GX - XX(j)
            BPY = GY - YY(j)
            BPZ = GZ - ZZ(j)
            BP = sqrt(BPX**2 + BPY**2 + BPZ**2)

            ! Calculate COSTH2 (cosine of angle between A and BP)
            COSTH2 = -(ABX * BPX + ABY * BPY + ABZ * BPZ) / (AB * BP)

            ! Calculate intermediate variables for velocity calculation
            V1 = APY * BPZ - APZ * BPY
            V2 = APZ * BPX - APX * BPZ
            V3 = APX * BPY - APY * BPX

            ! Calculate norm of area element
            H = sqrt(V1**2 + V2**2 + V3**2) / AB

            ! Skip calculation if any vector has zero magnitude or area element is zero
            if (AB <= 0.0 .or. AP <= 0.0 .or. BP <= 0.0 .or. H <= 0.0) then
                cycle
            end if

            ! Calculate magnitude of velocity
            VP = (COSTH1 + COSTH2) / (4.0 * pi * H)

            ! Calculate vector ABP (cross product of A and BP)
            ABPX = ABY * APZ - ABZ * APY
            ABPY = ABZ * APX - ABX * APZ
            ABPZ = ABX * APY - ABY * APX
            ABP = sqrt(ABPX**2 + ABPY**2 + ABPZ**2)

            ! Calculate components of velocity
            VPX=VP*ABPX/ABP; VPY=VP*ABPY/ABP; VPZ=VP*ABPZ/ABP

            VX = VX + VPX; VY = VY + VPY; VZ = VZ + VPZ
        end do
    end subroutine vortex

    subroutine velwak(gx,gy,gz,xvv,yvv,zvv)
        implicit none

        real(8), intent(in)     :: gx,gy,gz
        real(8), intent(inout)  :: xvv,yvv,zvv
        real(8)                 :: vvx,vvy,vvz
        real(8), dimension(4)   :: xx,yy,zz
        real(8) :: invz = 1.d0, invy= 1.d0


        integer :: NG,NA
        integer :: IK,JK

        NG = 0; invz = 1.d0 ! if NG =0 we do not invert the Z (for ground effect)
        xvv=0.d0; yvv=0.d0; zvv=0.d0

        do while (.true.) 
            do IK = 1, iteration
                do JK = 1, NPW
                    NA=0; invy = 1.d0 ! if NA = 0 we do not invert the Y (for symmetry)

                    xx = [x_wake(IK,ICW(JK,1)), x_wake(IK,ICW(JK,2)), x_wake(IK+1,ICW(JK,2)), x_wake(IK+1,ICW(JK,1))]
                    yy = [y_wake(IK,ICW(JK,1)), y_wake(IK,ICW(JK,2)), y_wake(IK+1,ICW(JK,2)), y_wake(IK+1,ICW(JK,1))]
                    zz = [z_wake(IK,ICW(JK,1)), z_wake(IK,ICW(JK,2)), z_wake(IK+1,ICW(JK,2)), z_wake(IK+1,ICW(JK,1))]

                    do while (.true.) !51
                        call vortex(gx,gy,gz,xx,yy,zz,vvx,vvy,vvz)
                        XVV=XVV+VVX*CW(IK,JK)*((-1)**NA)*(-1)**NG
                        YVV=YVV+VVY*CW(IK,JK)*((-1)**NA)*(-1)**NG
                        ZVV=ZVV+VVZ*CW(IK,JK)*((-1)**NA)*(-1)**NG
                        if ((config%idx_symm == 0).or.(NA == 1)) exit ! go to the next j panel
                            NA = 1; invy = - invy ! When NA is 1 we invert the Y (symmetry)
                    end do
                end do
            end do
            if ((config%idx_ground_eff == 0).or.(NG == 1)) exit
            NG = 1; invz = - invz ! When NG = 1 we invert the Z (for ground effect)
        end do
    end subroutine velwak

end program main