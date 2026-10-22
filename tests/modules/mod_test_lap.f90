module mod_test_lap
    use :: typedef, only:wp
    use :: initial_conditions, only:real_functor, compl_functor
    implicit none
    !width and wavevector small enough that O(h^2) FD error stays under tol
    real(wp), parameter :: alpha = 1.0e-4_wp
    real(wp), parameter :: kvec(3) = [2.0e-2_wp, 3.0e-2_wp, 1.5e-2_wp]

    type, extends(compl_functor) :: gaupsi0
    contains
        procedure :: eval => gauss_psi0
        procedure :: read_params => dummy_read
    end type gaupsi0

    type, extends(compl_functor) :: psi_analitic
    contains
        procedure :: eval => gauss_laplaceian
        procedure :: read_params => dummy_read2
    end type psi_analitic

    type(gaupsi0) :: dgpo
contains
    !initial conditions subroutine conforming to
    !abstract interface ./src/util/system
    !this is reasonably smooth and importantly has a known laplacian
    subroutine gauss_psi0(self, r, val)
        class(gaupsi0), intent(in) :: self
        real(wp), intent(in) :: r(3)
        complex(wp), intent(out) :: val

        if(.false.) then
            select type (s => self)
            end select
        end if

        val = exp(-alpha*sum(r**2))*exp(cmplx(0._wp, sum(kvec*r), kind=wp))
    end subroutine gauss_psi0

    !the analitical laplaceian of those initial conditions
    !abstract interface ./src/util/system
    subroutine gauss_laplaceian(self, r, val)
        class(psi_analitic), intent(in) :: self
        real(wp), intent(in) :: r(3)
        complex(wp), intent(out) :: val

        complex(wp) :: psi
        real(wp) :: r2, k2, kr

        if(.false.) then
            select type (s => self)
            end select
        end if

        r2 = sum(r**2)
        k2 = sum(kvec**2)
        kr = sum(kvec*r)
        call dgpo%eval(r, psi)
        ! nabla^2 [exp(-alpha r^2) exp(i k.r)]
        !   = [4 alpha^2 r^2 - 6 alpha - |k|^2 - 4 i alpha (k.r)] psi
        val = (4._wp*alpha**2*r2 - 6._wp*alpha - k2)*psi
        val = val - cmplx(0._wp, 4._wp*alpha*kr, kind=wp)*psi
    end subroutine gauss_laplaceian

    subroutine dummy_read(self, fname)
        class(gaupsi0), intent(inout) :: self
        character(len=*), intent(in) :: fname

        complex(wp) :: outval
        if(.false.) then
            write(*,*) fname
            call self%eval([1._wp, 1._wp, 1._wp], outval)
        end if
    end subroutine dummy_read
    subroutine dummy_read2(self, fname)
        class(psi_analitic), intent(inout) :: self
        character(len=*), intent(in) :: fname
        complex(wp) :: outval
        if(.false.) then
            write(*,*) fname
            call self%eval([1._wp, 1._wp, 1._wp], outval)
        end if
    end subroutine dummy_read2
end module mod_test_lap
