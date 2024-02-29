module settings
    implicit none
    type setting
        real(8) :: alpha = 0.d0, beta = 0.d0, gamma = 0.d0
        real(8) :: vinit = 0.d0, epsilon = 0.d0
        real(8) :: flight_height = 0.d0
        integer :: idx_symm = 0, idx_inch = 0, idx_ground_eff = 0
        real(8) :: dt = 0.03
        integer :: num_points = 0, num_panels = 0, num_triangles = 0, num_quads = 0
        integer :: mum_marked = 0
    end type setting
end module settings