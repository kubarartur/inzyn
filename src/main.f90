program main
    use :: typedef, only:wp
    use :: grid_util, only:grid_params
    use :: system_util, only:system
    use :: solver, only:H_GPe
    use :: io, only:write_arr
    use :: phys_const, only:pi => const_pi, m_to_a => conv_m_to_a0
    use :: initial_conditions, only:double_well, double_gaussian
    implicit none

    type(grid_params) :: gr
    type(system) :: sys
    type(H_GPe) :: H
    integer, parameter :: N(3) = [33, 35, 13]
    real(wp), parameter :: box_size(3) = [6e-6_wp, 6e-6_wp, 6e-6_wp]*m_to_a
    real(wp), parameter :: cell_size(3) = box_size/N

    call gr%init(N, cell_size)
    call sys%init(gr, double_gaussian)
    call H%init(gr, sys, double_well)

    call write_arr(sys%psi, gr, "Initial conditions for the wavefunction", "data/psi.txt")
    call write_arr(H%V_ext, gr, "The external trap potential", "data/V_ext.txt")

    write (*, *) "Hello, world!"
contains
end program main
