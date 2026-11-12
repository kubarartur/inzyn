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
    real(wp), allocatable :: r_right(:, :)
    real(wp), allocatable :: r_left(:, :)
    integer :: i
    integer :: halfx, endx, halfz
    integer :: argc
    character(len=:), allocatable :: out_dir_name
    integer :: stat

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
    allocate (r_right(3, 0:sol%real_steps_todo))
    allocate (r_left(3, 0:sol%real_steps_todo))

    endx = size(sys%psi, 1)
    halfx = size(sys%psi, 1)/2
    halfz = size(sys%psi, 3)/2

    call sol%find_E(gr, sys, H, energies(0))
    r_right(:, 0) = r_expct(gr, sys%psi(halfx:endx, :, halfz:halfz), [halfx, 0, halfz])
    r_left(:, 0) = r_expct(gr, sys%psi(1:halfx, :, halfz:halfz), [0, 0, halfz])

    do i = 1, sol%real_steps_todo
        write (*, *) "*****************************************************************************************"
        write (*, *) "real step: ", i, " / ", sol%real_steps_todo

        call sol%step_real(gr, sys, H, stat)
        if (stat /= 0) then
            call write_spacial(out_dir_name//"psi_final", sys%psi, gr%blc, gr%dr)
            call write_arr(out_dir_name//"E_n_of_n", energies)
            call write_arr(out_dir_name//"r_right_of_n", r_right)
            call write_arr(out_dir_name//"r_left_of_n", r_left)
            stop 1
        end if

        !cheap call, gives info as solver should hold E const
        call sol%find_E(gr, sys, H, energies(i), no_recalc=.true.)
        r_right(:, i) = r_expct(gr, sys%psi(halfx:endx, :, halfz:halfz), [halfx, 0, halfz])
        r_left(:, i) = r_expct(gr, sys%psi(1:halfx, :, halfz:halfz), [0, 0, halfz])

        if (modulo(i, sol%log_real_every) == 0) then
            block
                character(len=16) :: savename
                write (savename, '(A,I6.6)') "real_step_", i
                call write_spacial(out_dir_name//savename, sys%psi, gr%blc, gr%dr)
            end block
        end if
    end do

    call write_spacial(out_dir_name//"psi_final", sys%psi, gr%blc, gr%dr)
    call write_arr(out_dir_name//"E_n_of_n", energies)
    call write_arr(out_dir_name//"r_right_of_n", r_right)
    call write_arr(out_dir_name//"r_left_of_n", r_left)
contains

    function r_expct(gr_, psi, offset) result(ret)
        type(grid_params), intent(in) :: gr_
        complex(wp), intent(in) :: psi(:, :, :)
        integer, intent(in) :: offset(3)
        real(wp) :: ret(3)

        integer :: ii, jj, kk
        real(wp) :: loc_norm
        real(wp) :: loc_r(3)
        real(wp) :: dV
        real(wp) :: loc_n

        loc_norm = 0
        ret = 0
        dV = gr_%dr(1)*gr_%dr(2)*gr_%dr(3)

        do kk = 1, size(psi, 3)
        do jj = 1, size(psi, 2)
        do ii = 1, size(psi, 1)
            loc_r = gr_%blc + ([ii, jj, kk] + offset)*gr_%dr
            loc_n = dV*abs(psi(ii, jj, kk))**2
            loc_norm = loc_norm + loc_n

            ret = ret + loc_r*loc_n
        end do
        end do
        end do

        ret = ret/loc_norm
    end function r_expct
end program ground_seek
