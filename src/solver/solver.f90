module solver_mod
    use :: typedef, only:wp
    use :: grid_util, only:grid_params
    use :: system_util, only:system
    use :: hamiltonian, only:H_GPe
    use, intrinsic :: iso_fortran_env, only: stderr => error_unit

    implicit none
    private

    type, public :: solver
        real(wp) :: Euler_imag_dt
        !!timestep in forward Euler integrator for finding the gound state
        !!using Wick's rotation
        real(wp) :: time
        integer :: nsteps

    contains
        procedure :: init => init_solver
        procedure :: step_imag => imag_step_solver
    end type solver
contains

    subroutine init_solver(sol, dt)
        class(solver), intent(out) :: sol
        real(wp), intent(in) :: dt

        sol%nsteps = 0
        sol%time = 0._wp
        sol%Euler_imag_dt = dt
    end subroutine init_solver

    subroutine imag_step_solver(sol, gr, sys, H)
        class(solver), intent(inout) :: sol
        type(grid_params), intent(in) :: gr
        type(system), intent(inout) :: sys
        type(H_GPe), intent(inout) :: H

        complex(wp), allocatable :: psi_cp(:, :, :)
        real(wp) :: norm_curr

        !TODO: need for such call points to possible refactor
        !with the change the application of H is handled
        psi_cp = sys%psi

        !modify sys%psi in place
        call H%apply(gr, sys)

        !make the step in forward Euler
        sys%psi = psi_cp - sol%Euler_imag_dt*sys%psi
        sol%nsteps = sol%nsteps + 1
        sol%time = sol%time + sol%Euler_imag_dt

        !normalize the wavefunctoin
        call find_norm(sys%psi, gr%dr, norm_curr)

        !unsafe but if the norm is 0 then we don't have much to work with anyways
        sys%psi = sys%psi/norm_curr
    contains
        subroutine find_norm(arr, dr, norm)
            complex(wp), intent(in) :: arr(:, :, :)
            real(wp), intent(in) :: dr(3)
            real(wp), intent(out) :: norm

            integer :: ii, jj, kk
            real(wp) :: dV
            dV = dr(1)*dr(2)*dr(3)
            norm = 0

            do kk = 1, size(arr, 3)
            do jj = 1, size(arr, 2)
            do ii = 1, size(arr, 1)
                norm = norm + dV*arr(ii, jj, kk)*conjg(arr(ii, jj, kk))
            end do
            end do
            end do
            norm = sqrt(norm)
        end subroutine find_norm
    end subroutine imag_step_solver

end module solver_mod
