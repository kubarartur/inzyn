program main
    use :: typedef, only:wp
    use :: grid_util, only:grid
    use :: phys_const, only:pi => const_pi
    implicit none

    type(grid) :: gr
    real(wp) :: b(3)
    integer :: N(3)

    N= [12, 13, 15]
    b= [3._wp, pi, 22._wp/7._wp]
    call gr%init(N, b)

    write(*,*) "Hello, world!"
end program main
