module solver_mod
    use :: typedef, only:wp
    use :: grid_util, only:grid_params, norm
    use :: system_util, only:system
    use :: hamiltonian, only:H_GPe
    use, intrinsic :: iso_fortran_env, only: stderr => error_unit

    implicit none
    private

    type, public :: solver
        real(wp) :: Euler_imag_dt
        !!timestep in forward Euler integrator for finding the gound state
        !!using Wick's rotation
        real(wp) :: CN_real_dt
        !!real timestep in Crank-Nicolson scheme
        real(wp) :: time
        !!total time spent iterating
        integer :: nsteps
        !!current number of steps already done
        integer :: imag_steps_todo
        !! total number of steps to be done in the ground seek
        integer :: real_steps_todo
        !! total number of steps to be done in the real time evolution
        integer :: log_real_every
        !! number of steps to do before loggin a real time evolution step
    contains
        procedure :: init => init_solver
        procedure :: step_imag => imag_step_solver
        procedure :: step_real => real_step_solver
        procedure, nopass :: find_E => find_E_solver
    end type solver
contains

    subroutine init_solver(sol, idt, istd, rdt, rstd, lre)
        class(solver), intent(out) :: sol
        real(wp), intent(in) :: idt
        integer, intent(in) :: istd
        real(wp), intent(in) :: rdt
        integer, intent(in) :: rstd
        integer, intent(in) :: lre

        sol%nsteps = 0
        sol%time = 0._wp
        sol%Euler_imag_dt = idt
        sol%CN_real_dt = rdt
        sol%imag_steps_todo = istd
        sol%real_steps_todo = rstd
        sol%log_real_every = lre
    end subroutine init_solver

    subroutine imag_step_solver(sol, gr, sys, H, no_recalc, norm_first, delta_n)
        class(solver), intent(inout) :: sol
        type(grid_params), intent(in) :: gr
        type(system), intent(inout) :: sys
        type(H_GPe), intent(inout) :: H
        logical, optional, intent(in) :: no_recalc
        !! if true then don't recalculate the H%Vdd and H%n
        !! inside of H%apply
        logical, optional, intent(in) :: norm_first
        real(wp), optional, intent(out) :: delta_n
        !! if present track the maximum change of density

        complex(wp), allocatable :: psi_cp(:, :, :)

        if (present(norm_first)) then
            if (norm_first) then
                call norm(sys%psi, sys%ncz, gr%dr)
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

        !manual boundary conition enforcement once again
        sys%psi(:,:,1) = 0
        sys%psi(:,1,:) = 0
        sys%psi(1,:,:) = 0
        sys%psi(:,:,size(sys%psi, 3)) = 0
        sys%psi(:,size(sys%psi, 2),:) = 0
        sys%psi(size(sys%psi, 1),:,:) = 0

        call norm(sys%psi, sys%ncz, gr%dr)

        !after the step has been completed check how much the density
        !have changed in the cell it's changed the most
        if (present(delta_n)) then
            delta_n = maxval(abs(sys%psi)**2 - abs(psi_cp)**2)
        end if
    end subroutine imag_step_solver

    subroutine real_step_solver(sol, gr, sys, H)
        class(solver), intent(inout) :: sol
        type(grid_params), intent(in) :: gr
        type(system), intent(inout) :: sys
        type(H_GPe), intent(inout) :: H

        !set of predictors is needed to complete iteration
        type(system) :: sys_n
        type(system) :: sys_Hn
        type(system) :: sys_n1

        !helpers
        real(wp) :: dV, tol
        integer :: ii

        dV = gr%dr(1)*gr%dr(2)*gr%dr(3)
        tol = 1.e-8_wp

        associate (dt => sol%CN_real_dt)

            !this perseveres psi_n
            call sys_n%copy(sys)

            !this perseveres H_n\psi_n
            call sys_Hn%copy(sys)
            call H%apply(gr, sys_Hn)

            ii = 0
            do
                !this perseveres approximated \psi_{n+1}
                !needed only for convergence check
                call sys_n1%copy(sys)

                call H%apply(gr, sys)

                !atomic units \hbar omitted
                sys%psi = sys_n%psi - (0, 0.5_wp)*dt*(sys_Hn%psi + sys%psi)

                if (convergence_criterion(dV, sys%psi, sys_n1%psi, tol)) then
                    exit
                end if
                ii = ii + 1
            end do

        end associate

    contains

        function convergence_criterion(dV_, arr1, arr2, tol_) result(convergence)
            real(wp), intent(in) :: dV_
            complex(wp), intent(in) :: arr1(:, :, :)
            complex(wp), intent(in) :: arr2(:, :, :)
            real(wp), intent(in) :: tol_
            logical :: convergence

            integer :: i, j, k
            real(wp) :: diffn

            diffn = 0._wp
            do k = 1, size(arr1, 3)
            do j = 1, size(arr1, 2)
            do i = 1, size(arr1, 1)
                diffn = diffn + dV_*abs(arr1(i, j, k) - arr2(i, j, k))**2
            end do
            end do
            end do

            if (diffn < tol_**2) then
                convergence = .true.
            else
                convergence = .false.
            end if

            write (*, *) "-----------------------------------------"
            write (*, *) "CN point iter: ", ii
            write (*, *) "diffn: ", diffn
            write (*, *) "tol^2: ", tol**2
            write (*, *) "-----------------------------------------"

        end function convergence_criterion

    end subroutine real_step_solver

    subroutine find_E_solver(gr, sys, H, res, no_recalc)
        type(grid_params), intent(in) :: gr
        type(system), intent(in) :: sys
        type(H_GPe), intent(inout) :: H
        complex(wp), intent(out) :: res
        logical, optional, intent(in) :: no_recalc
            !! if true then don't recalculate the H%Vdd and H%n
            !! inside of H%apply

        if (present(no_recalc)) then
            res = H%find_E(gr, sys, no_recalc=no_recalc)
        else
            res = H%find_E(gr, sys)
        end if
    end subroutine find_E_solver

end module solver_mod
