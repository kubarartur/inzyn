module single_well_mod
    use :: typedef, only:wp
    use :: functor, only:real_functor
    use :: phys_const, only: &
        w => const_w, &
        m_to_a => conv_m_to_a0, &
        m => const_m, &
        nk_to_atm => conv_nK_to_hartree

    implicit none

    type, extends(real_functor) :: single_well
        real(wp) :: d = 1.5e-6_wp*m_to_a
        real(wp) :: alpha = 0.25_wp*m*(w(1)/1.5e-6_wp/m_to_a)**2
    contains
        procedure :: eval => eval_swell
        procedure :: read_params => read_swell_params
    end type

contains
    subroutine eval_swell(self, r, val)
        !!double well initial external trap potential
        !!see implementation for parameters
        class(single_well), intent(in) :: self
        real(wp), intent(in) :: r(3)
        real(wp), intent(out) :: val

        associate ( &
            d => self%d, &
            alpha => self%alpha)

            val = alpha*r(1)**4
            val = val + 0.5_wp*m*(w(2)*r(2))**2
            val = val + 0.5_wp*m*(w(3)*r(3))**2
        end associate
    end subroutine eval_swell

    subroutine read_swell_params(self, fname)
        class(single_well), intent(inout) :: self
        character(len=*), intent(in) :: fname

        integer :: io
        real(wp) :: d_temp

        real(wp) :: d = 1.5_wp*m_to_a
        namelist /single_well_parameters/ d

        d_temp = self%d

        io = 67
        open (io, file=trim(fname), status="old", action="read")
        read (io, nml=single_well_parameters)
        close (io)

        self%d = d*m_to_a*1e-6
        self%alpha = self%alpha*(d_temp/self%d)**2

        write (*, *) "Ex is: ", w(1)
    end subroutine read_swell_params
end module single_well_mod
