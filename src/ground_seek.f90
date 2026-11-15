program ground_seek
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

    real(wp), allocatable :: dns(:)
    integer :: i
    integer :: argc
    character(len=:), allocatable :: out_dir_name

    complex(wp), allocatable :: energies(:)
    complex(wp), allocatable :: contact(:)
    complex(wp), allocatable :: contact_ext(:)
    complex(wp), allocatable :: ext(:)
    complex(wp), allocatable :: kinetic(:)
    complex(wp), allocatable :: dipole(:)


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

    allocate (energies(0:sol%imag_steps_todo))
    allocate (contact(0:sol%imag_steps_todo))
    allocate (contact_ext(0:sol%imag_steps_todo))
    allocate (ext(0:sol%imag_steps_todo))
    allocate (kinetic(0:sol%imag_steps_todo))
    allocate (dipole(0:sol%imag_steps_todo))

    allocate (dns(sol%imag_steps_todo))

    call H%find_E(gr, sys, energies(0), &
                  contact=contact(0), contact_ext=contact_ext(0), &
                  ext=ext(0), kinetic=kinetic(0), dipole=dipole(0))

    do i = 1, sol%imag_steps_todo
        call sol%step_imag(gr, sys, H, no_recalc=.true., delta_n=dns(i))
        call H%find_E(gr, sys, energies(i), &
                      contact=contact(i), contact_ext=contact_ext(i), &
                      ext=ext(i), kinetic=kinetic(i), dipole=dipole(i))

        write (*, *) "Step: ", i, "/", sol%imag_steps_todo
        write (*, *) "Ediff: ", energies(i - 1) - energies(i)
    end do

    !"Result of the ground state seek
    call write_spacial(out_dir_name//"psi_final", sys%psi, gr%blc, gr%dr)
    !"Energy as function of imag time step"
    call write_arr(out_dir_name//"E_of_n", energies)
    call write_arr(out_dir_name//"contact_of_n", contact)
    call write_arr(out_dir_name//"contact_ext_of_n", contact_ext)
    call write_arr(out_dir_name//"ext_of_n", ext)
    call write_arr(out_dir_name//"kinetic_of_n", kinetic)
    call write_arr(out_dir_name//"dipole_of_n", dipole)

    !"Maximum change of the density in a single cell"
    call write_arr(out_dir_name//"delta_n_of_n", dns)
contains
end program ground_seek
