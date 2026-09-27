module solver
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
        procedure :: step => step_solver
    end type solver
contains

    subroutine init_solver(sol)
        class(solver), intent(out) :: sol

        sol%nsteps = 0
        sol%time= 0
    end subroutine init_solver

    subroutine step_solver(sol, gr, sys, H)
        class(solver), intent(in) :: sol
        type(grid_params), intent(in) :: gr
        type(system), intent(inout) :: sys
        type(H_GPe), intent(inout) :: H

    end subroutine step_solver
end module solver
