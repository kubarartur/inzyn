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

    complex(wp), allocatable :: energies(:)
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

    allocate (energies(0:sol%real_steps_todo))

    call sol%find_E(gr, sys, H, energies(0))

    do i = 1, sol%real_steps_todo
        write(*,*) "*****************************************************************************************"
        write(*,*) "real step: ", i , " / " , sol%real_steps_todo
        call sol%step_real(gr, sys, H)

        !cheap call, gives info as solver should hold E const
        call sol%find_E(gr, sys, H, energies(i), no_recalc=.true.)

        if( modulo(i, sol%log_real_every) == 0 ) then
            block
                character(len=16) :: savename
                write(savename, '(A,I6.6)') "real_step_", i
                call write_spacial(out_dir_name // savename, sys%psi, gr%blc, gr%dr)
            end block
        end if
    end do

    !"Result of the ground state seek
    call write_spacial(out_dir_name//"psi_final", sys%psi, gr%blc, gr%dr)
    call write_arr(out_dir_name//"E_n_of_n", energies)
contains
end program ground_seek
