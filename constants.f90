module constants
    implicit none
    real(8), parameter :: rho = 1.293           ! kg/m-3 T = 273K and P = 101.325 kPa
    real(8), parameter :: pi  = 4.d0*atan(1.d0) ! double precision PI
    real(8), parameter :: zero = 0.d0
    integer, parameter :: max_iterations = 60   ! The number of max iterations, and the number of max wakes

    character(len=*), dimension(3), parameter :: filenames = ['inputs.dat','points.dat','panels.dat'] ! Dont change the order
end module constants