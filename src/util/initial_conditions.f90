!> module conaining definitions of arbitrary functions
!! conforming to the to according interfaces
!! for using them to to populate grids within system
module initial_conditions
    use :: typedef, only:wp
    use :: phys_const, only:w => const_w, m => const_m, m_to_a => conv_m_to_a0
    implicit none

contains

    pure subroutine double_well(r, V)
        !!double well initial external trap potential
        !!see implementation for parameters
        real(wp), intent(in) :: r(3)
        real(wp), intent(out) :: V

        real(wp), parameter :: d = 1.5e-6_wp*m_to_a
        real(wp), parameter :: alpha = 0.25_wp*m*(w(1)/d)**2
        real(wp), parameter :: beta = 0.5_wp*m*w(1)**2

        V = alpha*r(1)**4 - beta*r(1)**2
        V = V + 0.5_wp*m*(w(2)*r(2))**2
        V = V + 0.5_wp*m*(w(3)*r(3))**2
    end subroutine double_well

    pure subroutine double_gaussian(r, phi)
        !!double gausian initial condition for wavefunction
        !!see implementation for parameters
        real(wp), intent(in) :: r(3)
        complex(wp), intent(out) :: phi

        real(wp), parameter :: d = 1.5e-6_wp*m_to_a
        real(wp), parameter :: sigma(3) = [0.5_wp*d, 0.5_wp*d, 0.9_wp*d]
        real(wp), parameter :: r01(3) = [-d, 0._wp, 0._wp]
        real(wp), parameter :: r02(3) = [d, 0._wp, 0._wp]
        !positions of the first and the second gaussian

        phi = exp(-0.5*sum(((r - r01)/sigma)**2)) + exp(-0.5*sum(((r - r02)/sigma)**2))
    end subroutine double_gaussian

    pure subroutine parser_debug_real(r, V)
        !!double well initial external trap potential
        !!see implementation for parameters
        real(wp), intent(in) :: r(3)
        real(wp), intent(out) :: V

       V = r(1) + r(2) + r(3)
    end subroutine parser_debug_real
end module initial_conditions
