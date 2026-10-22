!> module conaining definitions of arbitrary functions
!! conforming to the to according interfaces
!! for using them to to populate grids within system
module initial_conditions
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

    type, extends(real_functor), private :: double_well
        real(wp) :: d = 1.5e-6_wp*m_to_a
        real(wp) :: alpha = 0.25_wp*m*(w(1)/1.5e-6_wp/m_to_a)**2
        real(wp) :: beta = 0.5_wp*m*w(1)**2
    contains
        procedure :: eval => eval_dwell
        procedure :: read_params => read_dwell_params
    end type

    type, extends(compl_functor), private :: double_gaussian
        real(wp) :: d = 1.9e-6_wp*m_to_a
        real(wp) :: sigma(3) = [0.9_wp*1.9e-6_wp*m_to_a, 0.9_wp*1.9e-6_wp*m_to_a, 0.9_wp*1.9e-6_wp*m_to_a]
        real(wp) :: r01(3) = [-1.9e-6_wp*m_to_a, 0._wp, 0._wp]
        real(wp) :: r02(3) = [1.9e-6_wp*m_to_a, 0._wp, 0._wp]
    contains
        procedure :: eval => eval_dgauss
        procedure :: read_params => read_dgauss
    end type

    public :: get_real_functor
    public :: get_compl_functor
contains

    subroutine eval_dwell(self, r, val)
        !!double well initial external trap potential
        !!see implementation for parameters
        class(double_well), intent(in) :: self
        real(wp), intent(in) :: r(3)
        real(wp), intent(out) :: val

        associate ( &
            d => self%d, &
            alpha => self%alpha, &
            beta => self%beta &
            )

            val = alpha*r(1)**4 - beta*r(1)**2
            val = val + 0.5_wp*m*(w(2)*r(2))**2
            val = val + 0.5_wp*m*(w(3)*r(3))**2

        end associate
    end subroutine eval_dwell

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

    subroutine read_dwell_params(self, fname)
        class(double_well), intent(inout) :: self
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
        !self beta does not change upon different d
    end subroutine read_dwell_params

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

    function get_real_functor(funname) result(functor)
        character(len=*), intent(in) :: funname
        class(real_functor), allocatable :: functor

        if (trim(funname) == "double_well") then
            allocate (double_well :: functor)
        end if
    end function get_real_functor

    function get_compl_functor(funname) result(functor)
        character(len=*), intent(in) :: funname
        class(compl_functor), allocatable :: functor

        if (trim(funname) == "double_gaussian") then
            allocate (double_gaussian :: functor)
        end if
    end function get_compl_functor
end module initial_conditions
