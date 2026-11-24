module zeros_mod
    use :: typedef, only:wp
    use :: functor, only:real_functor
    use :: phys_const, only: &
        w => const_w, &
        m_to_a => conv_m_to_a0, &
        m => const_m, &
        nk_to_atm => conv_nK_to_hartree

    implicit none

    type, extends(real_functor) :: zeros
        integer :: no_unused = 0
    contains
        procedure :: eval => eval_zeros
        procedure :: read_params => read_zeros_params
    end type

contains
    subroutine eval_zeros(self, r, val)
        !!double well initial external trap potential
        !!see implementation for parameters
        class(zeros), intent(in) :: self
        real(wp), intent(in) :: r(3)
        real(wp), intent(out) :: val

        val = 0

        if(.false.) then
            write(*,*) r
            write(*,*) self%no_unused
        end if
    end subroutine eval_zeros

    subroutine read_zeros_params(self, fname)
        class(zeros), intent(inout) :: self
        character(len=*), intent(in) :: fname

        if (.false.) then
            block
                real(wp) :: dummy
                call self%eval([0._wp, 0._wp, 0._wp], dummy)
            end block
            write(*,*) fname
        end if
    end subroutine read_zeros_params
end module zeros_mod
