program test_dist
    use :: typedef, only:wp
    use :: grid_util, only:system
    use :: phys_const, only:m_to_a => conv_m_to_a0
    use :: initial_conditions, only:double_well, double_gaussian
    use, intrinsic :: iso_fortran_env, only: stderr => error_unit

    implicit none

    real(wp), parameter :: tol = 1e-11_wp

    type(system) :: gr
    integer, parameter :: N(3) = [33, 35, 13]
    real(wp), parameter :: box_size(3) = [6e-6_wp, 6e-6_wp, 6e-6_wp]*m_to_a
    real(wp), parameter :: cell_size(3) = box_size/N
    integer, parameter :: test_indices(3, 5) = reshape([1, 1, 1, &
                                                        2, 2, 2, &
                                                        5, 8, 3, &
                                                        11, 5, 2, &
                                                        2, 2, 2, &
                                                        2, 13, 7], [3, 5])

    call gr%init(N, cell_size, double_gaussian, double_well)

    block
        integer :: i, j
        integer :: a, b, c
        integer :: a1, a2, b1, b2, c1, c2
        real(wp) :: dist1, dist2
        real(wp) :: r1(3), r2(3)

        do i = 1, size(test_indices, 2)
        do j = i, size(test_indices, 2)
            r1 = gr%blc + gr%dr*test_indices(:, i)
            r2 = gr%blc + gr%dr*test_indices(:, j)
            dist1 = norm2(r1 - r2)

            ! test is not production code
            a1 = test_indices(1, i)
            a2 = test_indices(1, j)
            a = abs(a1 - a2)

            b1 = test_indices(2, i)
            b2 = test_indices(2, j)
            b = abs(b1 - b2)

            c1 = test_indices(3, i)
            c2 = test_indices(3, j)
            c = abs(c1 - c2)

            dist2 = gr%dist(a, b, c)
            if (abs(dist2 - dist1) > tol) then
                write (stderr, *) "for the pair of indices: "
                write (stderr, *) "i: ", test_indices(:, i)
                write (stderr, *) "j: ", test_indices(:, j)
                write (stderr, *) "The error of: ", abs(dist2 - dist1)
                write (stderr, *) "exceeded the tolerance of: ", tol
                stop 1
            end if
        end do
        end do

    end block
end program test_dist
