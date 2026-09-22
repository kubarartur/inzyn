program main
    use :: typedef, only:wp
    use :: grid_util, only:system
    use :: phys_const, only:pi => const_pi
    implicit none

    type(system) :: gr
    real(wp) :: b(3)
    integer :: N(3)

    N = [12, 13, 15]
    b = [3._wp, pi, 22._wp/7._wp]
    call gr%init(N, b, phi0, V_ext)

    write (*, *) "Hello, world!"
contains
    pure subroutine phi0(r, phi)
        real(wp), intent(in) :: r(3)
        complex(wp), intent(out) :: phi
        phi= 0._wp
    end subroutine phi0
    pure subroutine V_ext(r, phi)
        real(wp), intent(in) :: r(3)
        real(wp), intent(out) :: phi
        phi= 0._wp
    end subroutine V_ext
end program main
