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

    real(wp), allocatable :: r_right(:, :)
    real(wp), allocatable :: r_left(:, :)
    real(wp), allocatable :: r_right_max(:, :)
    real(wp), allocatable :: r_left_max(:, :)

    integer :: i
    integer :: halfx, endx, halfz
    integer :: argc
    character(len=:), allocatable :: out_dir_name
    integer :: stat

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

    allocate (energies(0:sol%real_steps_todo))
    allocate (contact(0:sol%real_steps_todo))
    allocate (contact_ext(0:sol%real_steps_todo))
    allocate (ext(0:sol%real_steps_todo))
    allocate (kinetic(0:sol%real_steps_todo))
    allocate (dipole(0:sol%real_steps_todo))

    allocate (r_right(3, 0:sol%real_steps_todo))
    allocate (r_left(3, 0:sol%real_steps_todo))
    allocate (r_right_max(3, 0:sol%real_steps_todo))
    allocate (r_left_max(3, 0:sol%real_steps_todo))

    endx = size(sys%psi, 1)
    halfx = size(sys%psi, 1)/2
    halfz = size(sys%psi, 3)/2

    call H%find_E(gr, sys, energies(0), &
                  contact=contact(0), contact_ext=contact_ext(0), &
                  ext=ext(0), kinetic=kinetic(0), dipole=dipole(0))


    call r_info(gr, sys%psi(halfx:endx, :, halfz:halfz), [halfx-1, 0, halfz],&
         r_expct=r_right(:, 0), r_max=r_right_max(:,0))


    call r_info(gr, sys%psi(1:halfx, :, halfz:halfz), [0, 0, halfz],&
        r_expct=r_left(:, 0), r_max=r_left_max(:,0))

    do i = 1, sol%real_steps_todo
        write (*, *) "*****************************************************************************************"
        write (*, *) "real step: ", i, " / ", sol%real_steps_todo

        call sol%step_real(gr, sys, H, stat)
        if (stat /= 0) then
            call write_spacial(out_dir_name//"psi_final", sys%psi, gr%blc, gr%dr)
            call write_arr(out_dir_name//"E_n_of_n", energies)

            call write_arr(out_dir_name//"r_right_of_n", r_right)
            call write_arr(out_dir_name//"r_left_of_n", r_left)

            call write_arr(out_dir_name//"r_right_max_of_n", r_right_max)
            call write_arr(out_dir_name//"r_left_max_of_n", r_left_max)
            stop 1
        end if

        call H%find_E(gr, sys, energies(i), &
                      contact=contact(i), contact_ext=contact_ext(i), &
                      ext=ext(i), kinetic=kinetic(i), dipole=dipole(i))


        call r_info(gr, sys%psi(halfx:endx, :, halfz:halfz), [halfx-1, 0, halfz],&
             r_expct=r_right(:, i), r_max=r_right_max(:,i))

        call r_info(gr, sys%psi(1:halfx, :, halfz:halfz), [0, 0, halfz],&
            r_expct=r_left(:, i), r_max=r_left_max(:,i))

        if (modulo(i, sol%log_real_every) == 0) then
            block
                character(len=16) :: savename
                write (savename, '(A,I6.6)') "real_step_", i
                call write_spacial(out_dir_name//savename, sys%psi, gr%blc, gr%dr)
            end block
        end if
    end do

    call write_spacial(out_dir_name//"psi_final", sys%psi, gr%blc, gr%dr)

    call write_arr(out_dir_name//"E_of_n", energies)
    call write_arr(out_dir_name//"contact_of_n", contact)
    call write_arr(out_dir_name//"contact_ext_of_n", contact_ext)
    call write_arr(out_dir_name//"ext_of_n", ext)
    call write_arr(out_dir_name//"kinetic_of_n", kinetic)
    call write_arr(out_dir_name//"dipole_of_n", dipole)

    call write_arr(out_dir_name//"r_right_of_n", r_right)
    call write_arr(out_dir_name//"r_left_of_n", r_left)

    call write_arr(out_dir_name//"r_right_max_of_n", r_right_max)
    call write_arr(out_dir_name//"r_left_max_of_n", r_left_max)
contains

    subroutine r_info(gr_, psi, offset, r_expct, r_max)
        type(grid_params), intent(in) :: gr_
        complex(wp), intent(in) :: psi(:, :, :)
        integer, intent(in) :: offset(3)
        real(wp), intent(out), optional :: r_expct(3)
        real(wp), intent(out), optional :: r_max(3)

        integer :: ii, jj, kk
        real(wp) :: loc_norm
        real(wp) :: dV

        real(wp) :: loc_r(3)
        real(wp) :: loc_n

        real(wp) :: max_n

        loc_norm = 0
        max_n = 0
        dV = gr_%dr(1)*gr_%dr(2)*gr_%dr(3)
        if(present(r_expct)) then
            r_expct = 0
        end if

        max_n = 0

        do kk = 1, size(psi, 3)
        do jj = 1, size(psi, 2)
        do ii = 1, size(psi, 1)

            loc_r = gr_%blc + ([ii, jj, kk] + offset)*gr_%dr
            loc_n = dV*abs(psi(ii, jj, kk))**2
            loc_norm = loc_norm + loc_n


            if (present(r_expct)) then
                r_expct = r_expct + loc_r*loc_n
            end if
            if (present(r_max)) then
                if (loc_n > max_n) then
                    r_max = loc_r
                    max_n = loc_n
                end if
            end if
        end do
        end do
        end do

        if (present(r_expct)) then
            r_expct = r_expct/loc_norm
        end if
    end subroutine r_info
end program ground_seek
