!> module conaining definitions of arbitrary functions
!! conforming to the to according interfaces
!! for using them to to populate grids within system
module functor
    use :: typedef, only:wp
    use :: phys_const, only:w => const_w, m => const_m, m_to_a => conv_m_to_a0
    implicit none

    type, abstract, public :: real_functor
    contains
        procedure(eval_real_interface), deferred, public :: eval
        procedure(read_real_params_interface), deferred, public :: read_params
    end type

    type, abstract, public :: compl_functor
    contains
        procedure(eval_compl_interface), deferred, public :: eval
        procedure(read_compl_params_interface), deferred, public :: read_params
    end type


    abstract interface
        subroutine eval_real_interface(self, r, val)
            use :: typedef
            import :: real_functor
            class(real_functor), intent(in) :: self
            real(wp), intent(in) :: r(3)
            real(wp), intent(out) :: val
        end subroutine eval_real_interface
        subroutine eval_compl_interface(self, r, val)
            use :: typedef
            import :: compl_functor
            class(compl_functor), intent(in) :: self
            real(wp), intent(in) :: r(3)
            complex(wp), intent(out) :: val
        end subroutine eval_compl_interface
        subroutine read_real_params_interface(self, fname)
            import :: real_functor
            class(real_functor), intent(inout) :: self
            character(len=*), intent(in) :: fname
        end subroutine
        subroutine read_compl_params_interface(self, fname)
            import :: compl_functor
            class(compl_functor), intent(inout) :: self
            character(len=*), intent(in) :: fname
        end subroutine
    end interface
end module functor
