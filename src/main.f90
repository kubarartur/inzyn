program main
    use :: typedef, only:wp
    use :: grid_util, only:grid_params, norm
    use :: system_util, only:system
    use :: hamiltonian, only:H_GPe
    use :: solver_mod, only:solver
    use :: initializer, only:innit
    use :: arr_iout, only:write_arr, write_spacial
    implicit none

    type(grid_params) :: gr
    type(system) :: sys
    type(H_GPe) :: H
    type(solver) :: sol

    complex(wp), allocatable :: energies(:)
    real(wp), allocatable :: dns(:)
    integer :: i
    integer :: argc
    character(len=:), allocatable :: out_dir_name

    argc = command_argument_count()
    if (argc < 1 .or. argc > 1) then
        write (*, *) "Usage: main <input/output_dir>"
        stop 1
    end if
    block
        character(len=255) :: temp
        call get_command_argument(1, temp)
        out_dir_name = trim(temp)
        call innit(trim(temp)//"input.nml", gr, sys, H, sol)
    end block

    call norm(sys%psi, sys%ncz, gr%dr)
    !"Initial conditions for the wavefunction"
    call write_spacial(out_dir_name//"psi_init", sys%psi, gr%blc, gr%dr)
    !"The external trap potential"
    call write_spacial(out_dir_name//"V_ext", H%V_ext, gr%blc, gr%dr)

    allocate (energies(0:sol%steps_todo))
    allocate (dns(sol%steps_todo))

    call sol%find_E(gr, sys, H, energies(0))

    do i = 1, sol%steps_todo
        call sol%step_imag(gr, sys, H, no_recalc=.true., delta_n=dns(i))
        call sol%find_E(gr, sys, H, energies(i))
        write (*, *) "Step: ", i, "/", sol%steps_todo
        write (*, *) "Ediff: ", energies(i - 1) - energies(i)
    end do

    !"Result of the ground state seek
    call write_spacial(out_dir_name//"psi_final", sys%psi, gr%blc, gr%dr)
    !"Energy as function of imag time step"
    call write_arr(out_dir_name//"E_of_n", energies)
    !"Maximum change of the density in a single cell"
    call write_arr(out_dir_name//"delta_n_of_n", dns)
contains
end program main
