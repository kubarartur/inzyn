!> module conaining definitions of arbitrary functions
!! conforming to the to according interfaces
!! for using them to to populate grids within system
module initial_conditions
    use :: typedef, only:wp
    use :: phys_const, only:w => const_w, m => const_m, m_to_a => conv_m_to_a0
    implicit none

    type, private :: double_well_params_
        real(wp) :: d = 1.5e-6_wp*m_to_a
        real(wp) :: alpha = 0.25_wp*m*(w(1)/1.5e-6_wp/m_to_a)**2
        real(wp) :: beta = 0.5_wp*m*w(1)**2
    contains
        procedure, public :: read_params => read_dwell_params
    end type

    type, private :: double_gaussian_params_
        real(wp) :: d = 1.9e-6_wp*m_to_a
        real(wp) :: sigma(3) = [0.9_wp*1.9e-6_wp*m_to_a, 0.9_wp*1.9e-6_wp*m_to_a, 0.9_wp*1.9e-6_wp*m_to_a]
        real(wp) :: r01(3) = [-1.9e-6_wp*m_to_a, 0._wp, 0._wp]
        real(wp) :: r02(3) = [1.9e-6_wp*m_to_a, 0._wp, 0._wp]
        !positions of the first and the second gaussian
    contains
        procedure, public :: read_params => read_dgaussian_params
    end type

    type(double_well_params_) :: double_well_params
    type(double_gaussian_params_) :: double_gaussian_params
contains

    pure subroutine double_well(r, V)
        !!double well initial external trap potential
        !!see implementation for parameters
        real(wp), intent(in) :: r(3)
        real(wp), intent(out) :: V

        associate ( &
            d => double_well_params%d, &
            alpha => double_well_params%alpha, &
            beta => double_well_params%beta &
            )

            V = alpha*r(1)**4 - beta*r(1)**2
            V = V + 0.5_wp*m*(w(2)*r(2))**2
            V = V + 0.5_wp*m*(w(3)*r(3))**2

        end associate
    end subroutine double_well

    pure subroutine double_gaussian(r, phi)
        !!double gausian initial condition for wavefunction
        !!see implementation for parameters
        real(wp), intent(in) :: r(3)
        complex(wp), intent(out) :: phi

        associate ( &
            d => double_gaussian_params%d, &
            sigma => double_gaussian_params%sigma, &
            r01 => double_gaussian_params%r01, &
            r02 => double_gaussian_params%r02 &
            )

            phi = exp(-0.5*sum(((r - r01)/sigma)**2)) + exp(-0.5*sum(((r - r02)/sigma)**2))
        end associate
    end subroutine double_gaussian

    pure subroutine parser_debug_real(r, V)
        !!double well initial external trap potential
        !!see implementation for parameters
        real(wp), intent(in) :: r(3)
        real(wp), intent(out) :: V

        V = r(1) + r(2) + r(3)
    end subroutine parser_debug_real

    subroutine read_dwell_params(self, fname)
        class(double_well_params_), intent(out) :: self
        character(len=*), intent(in) :: fname

        integer :: io
        real(wp) :: d_temp

        real(wp) :: d = 1.5_wp*m_to_a
        namelist /double_well_parameters/ d

        d_temp = self%d


        io = 67
        open (io, file=trim(fname), status="old", action="read")
        read (io, nml=double_well_parameters)
        close (io)

        self%d = d*m_to_a*1e-6
        self%alpha = self%alpha*d_temp/self%d
    end subroutine read_dwell_params

    subroutine read_dgaussian_params(self, fname)
        class(double_gaussian_params_), intent(out) :: self
        character(len=*), intent(in) :: fname

        integer :: io

        namelist /double_gaussian_parameters/ d, sigma, r01, r02

        real(wp) :: d = 1.9_wp*m_to_a
        real(wp) :: sigma(3) = [0.9_wp*1.9_wp*m_to_a, 0.9_wp*1.9_wp*m_to_a, 0.9_wp*1.9_wp*m_to_a]
        real(wp) :: r01(3) = [-1.9_wp*m_to_a, 0._wp, 0._wp]
        real(wp) :: r02(3) = [1.9_wp*m_to_a, 0._wp, 0._wp]

        io = 67
        open (io, file=trim(fname), status="old", action="read")
        read (io, nml=double_gaussian_parameters)
        close (io)

        self%d = d*m_to_a*1e-6
        self%sigma = sigma*m_to_a*1e-6
        self%r01 = r01*m_to_a*1e-6
        self%r02 = r02*m_to_a*1e-6
    end subroutine read_dgaussian_params
end module initial_conditions
