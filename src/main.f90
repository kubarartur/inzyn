program main
    use :: typedef, only:wp
    use :: grid_util, only:grid_params
    use :: system_util, only:system
    use :: hamiltonian, only:H_GPe
    use :: solver_mod, only:solver
    use :: io, only:write_arr
    use :: phys_const, only:m_to_a => conv_m_to_a0
    use :: initial_conditions, only:double_well, double_gaussian
    implicit none

    type(grid_params) :: gr
    type(system) :: sys
    type(H_GPe) :: H
    type(solver) :: sol

    integer :: N(3) = [50, 50, 40]
    !!the size of the simulation grid
    integer :: Nstep = 1000
    !!number of steps of the solver in Wick's rotation part
    real(wp) :: box_size(3) = [6e-6_wp, 6e-6_wp, 6e-6_wp]
    !!size of the simulation box. Here in meters for consistency with namelist.
    !!unit conversion later on
    real(wp) :: dtau = 1e-5_wp
    !!timestep in the Wick's rotation

    real(wp) :: cell_size(3)
    complex(wp), allocatable :: energies(:)
    integer :: i
    integer :: io
    character(len=:), allocatable :: out_dir_name
    integer :: argc
    namelist /solver_params/ N, Nstep, box_size, dtau

    !read the sovler params from file
    io = 67
    open (io, file="input.nml", status="old", action="read")
    read (io, nml=solver_params)
    close (io)

    argc = command_argument_count()
    if (argc < 1 .or. argc > 1) then
        write (*, *) "Usage: main <output_dir>"
        stop 1
    end if
    block
        character(len=255) :: temp
        call get_command_argument(1, temp)
        out_dir_name = trim(temp)
    end block

    box_size = box_size*m_to_a
    cell_size = box_size/N

    call gr%init(N, cell_size)
    call sys%init(gr, double_gaussian)
    call H%init(gr, sys, double_well)
    call sol%init(dtau)

    call write_arr(sys%psi, gr, "Initial conditions for the wavefunction", out_dir_name//"psi_init.txt")
    call write_arr(H%V_ext, gr, "The external trap potential", out_dir_name//"V_ext.txt")

    allocate (energies(Nstep))
    do i = 1, Nstep
        call sol%step_imag(gr, sys, H)
        call sol%find_E(gr, sys, H, energies(i))
        write (*, *) "Step: ", i, "/", Nstep
    end do

    call write_arr(sys%psi, gr, "Result of the ground state seek", out_dir_name//"psi_final.txt")
    call write_arr(energies, [0._wp, 1._wp*Nstep], "Energy as function of imag time step", &
                   out_dir_name//"E_of_n.txt")
contains
end program main
