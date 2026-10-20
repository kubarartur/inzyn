program main
    use :: typedef, only:wp
    use :: grid_util, only:grid_params, norm
    use :: system_util, only:system
    use :: hamiltonian, only:H_GPe
    use :: solver_mod, only:solver
    use :: io, only:write_arr
    use :: phys_const, only:m_to_a => conv_m_to_a0
    use :: initial_conditions, only: &
        double_well, &
        double_gaussian, &
        double_well_params, &
        double_gaussian_params

    implicit none

    type(grid_params) :: gr
    type(system) :: sys
    type(H_GPe) :: H
    type(solver) :: sol

    integer :: N(3) = [50, 50, 40]
    !!the size of the simulation grid
    integer :: Nstep = 100
    !!number of steps of the solver in Wick's rotation part
    real(wp) :: box_size(3) = [6e-6_wp, 6e-6_wp, 6e-6_wp]
    !!size of the simulation box. Here in meters for consistency with namelist.
    !!unit conversion later on
    real(wp) :: dtau = 1e-5_wp
    !!timestep in the Wick's rotation
    integer :: ncz = 10000
    !!number of particles in the system

    real(wp) :: cell_size(3)
    complex(wp), allocatable :: energies(:)
    integer :: i
    integer :: iout
    character(len=:), allocatable :: out_dir_name
    integer :: argc
    namelist /solver_params/ N, Nstep, box_size, dtau, ncz

    argc = command_argument_count()
    if (argc < 1 .or. argc > 1) then
        write (*, *) "Usage: main <input/output_dir>"
        stop 1
    end if
    block
        character(len=255) :: temp
        call get_command_argument(1, temp)
        out_dir_name = trim(temp)
        !read the sovler params from file
        iout = 67
        open (iout, file=trim(temp)//"input.nml", status="old", action="read")
        read (iout, nml=solver_params)
        close (iout)
        call double_gaussian_params%read_params(trim(temp)//"input.nml")
        call double_well_params%read_params(trim(temp)//"input.nml")
    end block

    box_size = box_size*m_to_a
    cell_size = box_size/N

    call gr%init(N, cell_size)
    call sys%init(gr, double_gaussian, ncz)
    call H%init(gr, sys, double_well)
    call sol%init(dtau)

    call norm(sys%psi, sys%nczr, gr%dr)
    call write_arr(sys%psi, gr, "Initial conditions for the wavefunction", out_dir_name//"psi_init.txt")
    call write_arr(H%V_ext, gr, "The external trap potential", out_dir_name//"V_ext.txt")

    allocate (energies(0:Nstep))

    call sol%find_E(gr, sys, H, energies(0))

    do i = 1, Nstep
        call sol%step_imag(gr, sys, H, no_recalc=.true.)
        call sol%find_E(gr, sys, H, energies(i))
        write (*, *) "Step: ", i, "/", Nstep
        write (*, *) "Ediff: ", energies(i-1) - energies(i)
    end do

    call write_arr(sys%psi, gr, "Result of the ground state seek", out_dir_name//"psi_final.txt")
    call write_arr(energies, [0._wp, 1._wp*Nstep], "Energy as function of imag time step", &
                   out_dir_name//"E_of_n.txt")
contains
end program main
