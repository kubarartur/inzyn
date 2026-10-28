module initializer
    use :: typedef, only:wp
    use :: phys_const, only:m_to_a => conv_m_to_a0
    use :: grid_util, only:grid_params
    use :: solver_mod, only:solver
    use :: system_util, only:system
    use :: hamiltonian, only:H_GPe
    use :: functor, only:real_functor, compl_functor
    use :: double_gaussian_mod, only:double_gaussian
    use :: double_well_mod, only:double_well
    use :: from_file_mod, only:from_file
    implicit none

contains

    subroutine innit(fname, gr, sys, H, sol)
        character(len=*), intent(in) :: fname
        type(grid_params), intent(out) :: gr
        type(system), intent(out) :: sys
        type(H_GPe), intent(out) :: H
        type(solver), intent(out) :: sol

        call init_grid()
        call init_sys()
        call init_H()
        call init_solver()
    contains

        subroutine init_grid()
            integer :: iout
            integer :: N(3) = [50, 50, 30]
            real(wp) :: box_size(3) = [5., 5., 5.] !in micrometers
            real(wp) :: bottom_left_corner(3) = [-0.67_wp, -0.67_wp, -0.67_wp] !warning special value
            real(wp) :: dr(3)
            namelist /grid_params/ N, box_size, bottom_left_corner

            iout = 67
            open (iout, file=trim(fname), status="old", action="read")
            read (iout, nml=grid_params)
            close (iout)

            box_size = box_size*m_to_a*1e-6
            dr = box_size/N

            if (maxval(abs(bottom_left_corner - [-0.67_wp, -0.67_wp, -0.67_wp])) < 1e-20) then
                call gr%init(N, dr)
            else
                call gr%init(N, dr, bottom_left_corner=bottom_left_corner)
            end if
        end subroutine init_grid

        subroutine init_sys()
            integer :: iout
            class(compl_functor), allocatable :: psi0

            integer :: n_part = 20000
            logical :: load_sys = .false.
            character(len=255) :: load_fname = ""
            character(len=255) :: psi0_type = ""
            namelist /system_params/ n_part, load_sys, load_fname, psi0_type

            iout = 67
            open (iout, file=trim(fname), status="old", action="read")
            read (iout, nml=system_params)
            close (iout)

            if (load_sys) then
                psi0 = get_compl_functor("preloaded")
                call psi0%read_params(trim(fname))
                call sys%init(gr, psi0, n_part)
            else
                psi0 = get_compl_functor(trim(psi0_type))
                call psi0%read_params(trim(fname))
                call sys%init(gr, psi0, n_part)
            end if
        end subroutine init_sys

        subroutine init_solver()
            integer :: iout
            real(wp) :: Euler_imag_dt = 1e8
            integer :: steps_todo = 3
            namelist /solver_params/ Euler_imag_dt, steps_todo

            iout = 67
            open (iout, file=trim(fname), status="old", action="read")
            read (iout, nml=solver_params)
            close (iout)

            call sol%init(Euler_imag_dt, steps_todo)
        end subroutine init_solver

        subroutine init_H()
            integer :: iout
            class(real_functor), allocatable :: V_ex

            character(len=255) :: V_ext_type
            namelist /H_params/ V_ext_type

            iout = 67
            open (iout, file=trim(fname), status="old", action="read")
            read (iout, nml=H_params)
            close (iout)

            V_ex = get_real_functor(trim(V_ext_type))
            call V_ex%read_params(trim(fname))
            call H%init(gr, sys, V_ex)
        end subroutine init_H
    end subroutine innit

    function get_real_functor(funname) result(func)
        character(len=*), intent(in) :: funname
        class(real_functor), allocatable :: func

        if (trim(funname) == "double_well") then
            allocate (double_well :: func)
        end if
    end function get_real_functor

    function get_compl_functor(funname) result(func)
        character(len=*), intent(in) :: funname
        class(compl_functor), allocatable :: func

        if (trim(funname) == "double_gaussian") then
            allocate (double_gaussian :: func)
        else if (trim(funname) == "preloaded") then
            allocate (from_file :: func)
        end if
    end function get_compl_functor
end module initializer
