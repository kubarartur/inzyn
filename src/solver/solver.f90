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
        !!total time spent iterating
        integer :: nsteps
        !!current number of steps already done
        integer :: steps_todo
        !! total number of steps to be done

    contains
        procedure :: init => init_solver
        procedure :: step_imag => imag_step_solver
        procedure, nopass :: find_E => find_E_solver
    end type solver
contains

    subroutine init_solver(sol, dt, std)
        class(solver), intent(out) :: sol
        real(wp), intent(in) :: dt
        integer, intent(in) :: std

        sol%nsteps = 0
        sol%time = 0._wp
        sol%Euler_imag_dt = dt
        sol%steps_todo = std
    end subroutine init_solver

    subroutine imag_step_solver(sol, gr, sys, H, no_recalc, norm_first)
        class(solver), intent(inout) :: sol
        type(grid_params), intent(in) :: gr
        type(system), intent(inout) :: sys
        type(H_GPe), intent(inout) :: H
        logical, optional, intent(in) :: no_recalc
        logical, optional, intent(in) :: norm_first
        !! if true then don't recalculate the H%Vdd and H%n
        !! inside of H%apply

        complex(wp), allocatable :: psi_cp(:, :, :)
        real(wp) :: norm_curr

        if (present(norm_first)) then
            if (norm_first) then
                call find_norm(sys%psi, gr%dr, norm_curr)
                sys%psi = sqrt(sys%nczr)*sys%psi/norm_curr
            end if
        end if
        !TODO: need for such call points to possible refactor
        !with the change the application of H is handled
        allocate (psi_cp( &
                  size(sys%psi, 1), &
                  size(sys%psi, 2), &
                  size(sys%psi, 3) &
                  ))
        psi_cp = sys%psi

        if (present(no_recalc)) then
            call H%apply(gr, sys, no_recalc=no_recalc)
        else
            call H%apply(gr, sys)
        end if

        sys%psi = psi_cp - sol%Euler_imag_dt*sys%psi
        sol%nsteps = sol%nsteps + 1
        sol%time = sol%time + sol%Euler_imag_dt

        !normalize the wavefunctoin
        call find_norm(sys%psi, gr%dr, norm_curr)

        sys%psi = sqrt(sys%nczr)*sys%psi/norm_curr
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
                norm = norm + dV*abs(arr(ii, jj, kk))**2
            end do
            end do
            end do
            norm = sqrt(norm)
        end subroutine find_norm
    end subroutine imag_step_solver

    subroutine find_E_solver(gr, sys, H, res, no_recalc)
        type(grid_params), intent(in) :: gr
        type(system), intent(in) :: sys
        type(H_GPe), intent(inout) :: H
        complex(wp), intent(out) :: res
        logical, optional, intent(in) :: no_recalc
            !! if true then don't recalculate the H%Vdd and H%n
            !! inside of H%apply

        if (present(no_recalc)) then
            res= H%find_E(gr, sys, no_recalc=no_recalc)
        else
            res= H%find_E(gr, sys)
        end if
    end subroutine find_E_solver

end module solver_mod
