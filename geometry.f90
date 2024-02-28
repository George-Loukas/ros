module geometry
    use constants
    implicit none
    type panel_
        real(8), dimension(4,3) :: point = 0.d0
        integer, dimension(4)   :: mark = 0
        real(8), dimension(3)   :: centroid
        real(8), dimension(3)   :: normal
        real(8), dimension(3)   :: long
        real(8), dimension(3)   :: tang
        real(8)                 :: ds
        integer                 :: section = 0
        logical                 :: is_triangle = .false.

    contains
        procedure :: move_panel       => move_panel
        procedure :: get_centroid     => get_centroid
        procedure :: get_vectors      => get_vectors
    end type panel_
    contains


    subroutine move_panel(this,yaw, pitch, roll,height_flight,idx_ground_eff)
        use settings
        use constants
        implicit none
        class(panel_), intent(inout) :: this
        !! LOCAL VARIABLES
        integer :: i
        integer :: idx_ground_eff
        real(8) :: yaw, pitch, roll
        real(8) :: height_flight
        real(8) :: yawMat(3,3), pitchMat(3,3), rollMat(3,3), rotMat(3,3)
        real(8) :: sina,cosa,sinb,cosb,sing,cosg, temp_x, temp_y, temp_z

        yawMat  = 0; pitchMat = 0; rollMat = 0

        sina = sin(yaw)
        cosa = cos(yaw)

        sinb = sin(pitch)
        cosb = cos(pitch)

        sing = sin(roll)
        cosg = cos(roll)

        yawMat     = reshape([ cosa, -sina, 0.d0, &
                               sina,  cosa, 0.d0, &
                               0.d0,  0.d0, 1.d0], shape(yawMat))

        pitchMat   = reshape([ cosb, 0.d0, sinb, &
                               0.d0, 1.d0, 0.d0, &
                              -sinb, 0.d0, cosb], shape(pitchMat))

        rollMat    = reshape([ 1.d0, 0.d0,   0.d0,  &
                               0.d0, cosg, -sing, &
                               0.d0, sing,  cosg], shape(yawMat))

        rotMat     = matmul(matmul(yawMat,pitchMat),rollMat)


        do i = 1,4
            temp_x = this%point(i,1)
            temp_y = this%point(i,2)
            temp_z = this%point(i,3)

            this%point(i,1) = (rotMat(1,1)*temp_x + rotMat(1,2)*temp_y + rotMat(1,3)*temp_z) + height_flight*idx_ground_eff
            this%point(i,2) = (rotMat(2,1)*temp_x + rotMat(2,2)*temp_y + rotMat(2,3)*temp_z) + height_flight*idx_ground_eff
            this%point(i,3) = (rotMat(3,1)*temp_x + rotMat(3,2)*temp_y + rotMat(3,3)*temp_z) + height_flight*idx_ground_eff
        end do
    end subroutine move_panel

    subroutine get_centroid(this)
        implicit none
        class(panel_), intent(inout) :: this
        if (this%is_triangle) then
            this%centroid(1) = (this%point(1,1) + this%point(2,1) + this%point(3,1))/3.d0
            this%centroid(2) = (this%point(1,2) + this%point(2,2) + this%point(3,2))/3.d0
            this%centroid(3) = (this%point(1,3) + this%point(2,3) + this%point(3,3))/3.d0
        else
            this%centroid(1) = (this%point(1,1) + this%point(2,1) + this%point(3,1) + this%point(4,1))/4.d0
            this%centroid(2) = (this%point(1,2) + this%point(2,2) + this%point(3,2) + this%point(4,2))/4.d0
            this%centroid(3) = (this%point(1,3) + this%point(2,3) + this%point(3,3) + this%point(4,3))/4.d0
        end if
    end subroutine get_centroid

    subroutine get_vectors(this)
        implicit none
        class(panel_) :: this
        real(8) :: x13,y13,z13, x42,y42,z42, xx, yy, zz
        real(8) :: x12,y12,z12, x14,y14,z14
        real(8) :: xm,ym,zm
        real(8) :: norm, ds1, ds2

        x13 = this%point(1,1) - this%point(3,1)
        y13 = this%point(1,2) - this%point(3,2)
        z13 = this%point(1,3) - this%point(3,3)

        x42 = this%point(4,1) - this%point(2,1)
        y42 = this%point(4,2) - this%point(2,2)
        z42 = this%point(4,3) - this%point(2,3)

        xx = y13*z42 - z13*y42
        yy = z13*x42 - x13*z42
        zz = x13*y42 - y13*z42

        norm = sqrt(xx**2 + yy**2 + zz**2)

        this%normal(1) = - xx/norm
        this%normal(2) = - yy/norm
        this%normal(3) = - zz/norm

        x12 = this%point(2,1) - this%point(1,1)
        y12 = this%point(2,2) - this%point(1,2)
        z12 = this%point(2,3) - this%point(1,3)

        x14 = this%point(4,1) - this%point(1,1)
        y14 = this%point(4,2) - this%point(1,2)
        z14 = this%point(4,3) - this%point(1,3)

        xx = y13*z12 - z13*y12
        yy = z13*x12 - x13*z12
        zz = x13*y12 - y13*x12

        ds1 = 0.5*sqrt(xx**2 + yy**2 + zz**2)

        xx = y13*z14 - z13*y14
        yy = z13*x14 - x13*z14
        zz = x13*y14 - y13*x14

        ds2 = 0.5*sqrt(xx**2 + yy**2 + zz**2)

        this%ds = ds1 + ds2

        xm = (this%point(1,1) + this%point(2,1))/2.d0
        ym = (this%point(1,2) + this%point(2,2))/2.d0
        zm = (this%point(1,3) + this%point(2,3))/2.d0

        xx = xm - this%centroid(1)
        yy = ym - this%centroid(2)
        zz = zm - this%centroid(3)
        
        norm = sqrt(xx**2 + yy**2 + zz**2)

        this%long(1) = xx/norm
        this%long(2) = yy/norm
        this%long(3) = zz/norm

        xx = this%long(2)*this%normal(3) - this%long(3)*this%normal(2) 
        yy = this%long(3)*this%normal(1) - this%long(1)*this%normal(3)
        zz = this%long(1)*this%normal(2) - this%long(2)*this%normal(1)

        norm = sqrt(xx**2 + yy**2 + zz**2)

        this%tang(1) = xx/norm
        this%tang(2) = yy/norm
        this%tang(3) = zz/norm
    end subroutine get_vectors
end module geometry