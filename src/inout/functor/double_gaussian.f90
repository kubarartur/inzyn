module double_gaussian_mod
    use :: typedef, only:wp
    use :: functor, only:compl_functor
    use :: phys_const, only:m_to_a => conv_m_to_a0
    implicit none

    type, extends(compl_functor) :: double_gaussian
        real(wp) :: d = 1.9e-6_wp*m_to_a
        real(wp) :: sigma(3) = [0.9_wp*1.9e-6_wp*m_to_a, 0.9_wp*1.9e-6_wp*m_to_a, 0.9_wp*1.9e-6_wp*m_to_a]
        real(wp) :: r01(3) = [-1.9e-6_wp*m_to_a, 0._wp, 0._wp]
        real(wp) :: r02(3) = [1.9e-6_wp*m_to_a, 0._wp, 0._wp]
    contains
        procedure :: eval => eval_dgauss
        procedure :: read_params => read_dgauss
    end type

contains
    subroutine eval_dgauss(self, r, val)
        !!double gausian initial condition for wavefunction
        !!see implementation for parameters
        class(double_gaussian), intent(in) :: self
        real(wp), intent(in) :: r(3)
        complex(wp), intent(out) :: val

        associate ( &
            d => self%d, &
            sigma => self%sigma, &
            r01 => self%r01, &
            r02 => self%r02 &
            )

            val = exp(-0.5*sum(((r - r01)/sigma)**2)) + exp(-0.5*sum(((r - r02)/sigma)**2))
        end associate
    end subroutine eval_dgauss

    subroutine read_dgauss(self, fname)
        class(double_gaussian), intent(inout) :: self
        character(len=*), intent(in) :: fname

        integer :: io

        real(wp) :: d = 1.9_wp*m_to_a*1e-6
        real(wp) :: sigma(3) = &
                    [0.9_wp*1.9_wp*m_to_a*1e-6, 0.9_wp*1.9_wp*m_to_a*1e-6, 0.9_wp*1.9_wp*m_to_a*1e-6]
        real(wp) :: r01(3) = [-1.9_wp*m_to_a*1e-6, 0._wp, 0._wp]
        real(wp) :: r02(3) = [1.9_wp*m_to_a*1e-6, 0._wp, 0._wp]

        namelist /double_gaussian_parameters/ d, sigma, r01, r02

        io = 67
        open (io, file=trim(fname), status="old", action="read")
        read (io, nml=double_gaussian_parameters)
        close (io)

        self%d = d*m_to_a*1e-6
        self%sigma = sigma*m_to_a*1e-6
        self%r01 = r01*m_to_a*1e-6
        self%r02 = r02*m_to_a*1e-6
    end subroutine read_dgauss
end module double_gaussian_mod
