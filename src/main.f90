program main
    use :: typedef, only:wp
    use :: grid_util, only:system
    use :: phys_const, only:pi => const_pi , m_to_a => conv_m_to_a0
    use :: initial_conditions, only:double_well, double_gaussian
    implicit none

    type(system) :: gr
    integer, parameter :: N(3) = [13, 15, 7]
    real(wp), parameter :: box_size(3) = [6e-6_wp, 6e-6_wp, 6e-6_wp]*m_to_a
    real(wp), parameter :: cell_size(3) = box_size/N

    call gr%init(N, cell_size, double_gaussian, double_well)
    call gr%dump_psi("data/psi.txt")
    call gr%dump_v_ext("data/V_ext.txt")

    write (*, *) "Hello, world!"
contains
end program main
