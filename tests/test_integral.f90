program test_integral
    use :: typedef, only:wp
    use :: grid_util, only:grid_params
    use :: hamiltonian, only:vdd_boundary, laplaceian
    use :: phys_const, only:pi => const_pi
    use, intrinsic :: iso_fortran_env, only: stderr => error_unit
    implicit none

    type(grid_params) :: gr

    real(wp), parameter :: cell_size(3) = [2._wp, 0.7_wp, 3._wp]
    integer, parameter :: N(3) = [60, 60, 60]
    real(wp), parameter :: tol = 1._wp
    real(wp) :: max_err

    real(wp), allocatable :: source_arr(:, :, :)
    real(wp), allocatable :: res_num(:, :, :)
    real(wp), allocatable :: res_analitic(:, :, :)
    real(wp) :: r0(3)
    real(wp) :: sigma

    integer :: i, j, k

    call gr%init(n, cell_size)

    r0 = gr%blc + N*gr%dr*[0.45_wp, 0.5_wp, 0.6_wp]
    sigma = 3._wp*gr%dr(1)

    allocate (source_arr(N(1), N(2), N(3)))
    allocate (res_num(N(1), N(2), N(3)))
    allocate (res_analitic(N(1), N(2), N(3)))


    do k = 1, N(1)
    do j = 1, N(2)
    do i = 1, N(3)
        call gaussian(gr%blc + [k, j, i]*gr%dr, source_arr(k, j, i))
    end do
    end do
    end do

    do i = 1, N(3), N(3) - 1
    do j = 1, gr%dims(2)
    do k = 1, gr%dims(1)
        call res_errf(gr%blc + [k, j, i]*gr%dr, res_analitic(k, j, i))
    end do
    end do
    end do


    do i = 1, gr%dims(3)
    do j = 1, N(2), N(2) - 1
    do k = 1, gr%dims(1)
        call res_errf(gr%blc + [k, j, i]*gr%dr, res_analitic(k, j, i))
    end do
    end do
    end do


    do i = 1, gr%dims(3)
    do j = 1, gr%dims(2)
    do k = 1, N(1), N(1) - 1
        call res_errf(gr%blc + [k, j, i]*gr%dr, res_analitic(k, j, i))
    end do
    end do
    end do

    call vdd_boundary(res_num, gr, source_arr)

    block
        integer :: imax, jmax, kmax
        real(wp) :: err

        max_err = 0._wp
        imax = 1
        jmax = 1
        kmax = 1
        do i = 1, N(3), N(3) - 1
        do j = 1, gr%dims(2)
        do k = 1, gr%dims(1)

            err = abs(res_num(k, j, i) - res_analitic(k, j, i))/abs(res_analitic(k, j, i) + 1e-8)
            if (err > max_err) then
                max_err = err
                kmax = k
                jmax = j
                imax = i
            end if

        end do
        end do
        end do

        do i = 1, gr%dims(3)
        do j = 1, N(2), N(2) - 1
        do k = 1, gr%dims(1)

            err = abs(res_num(k, j, i) - res_analitic(k, j, i))/abs(res_analitic(k, j, i) + 1e-8)
            if (err > max_err) then
                max_err = err
                kmax = k
                jmax = j
                imax = i
            end if

        end do
        end do
        end do

        do i = 1, gr%dims(3)
        do j = 1, gr%dims(2)
        do k = 1, N(1), N(1) - 1

            err = abs(res_num(k, j, i) - res_analitic(k, j, i))/abs(res_analitic(k, j, i) + 1e-8)
            if (err > max_err) then
                max_err = err
                kmax = k
                jmax = j
                imax = i
            end if

        end do
        end do
        end do
        if (max_err > tol) then
            write (stderr, *) "integral mismatch at index: ", kmax, jmax, imax
            write (stderr, *) "numerical:  ", res_num(kmax, jmax, imax)
            write (stderr, *) "analitic:  ", res_analitic(kmax, jmax, imax)
            write (stderr, *) "The error of: ", max_err
            write (stderr, *) "exceeded the tolerance of: ", tol
            stop 1
        end if
    end block
    write (*, *) "Passed with max error of: ", max_err, "for tol of: ", tol

contains

    pure subroutine gaussian(r, val)
        real(wp), intent(in) :: r(3)
        real(wp), intent(out) :: val

        real(wp) :: d
        real(wp) :: sigma_
        real(wp) :: r0_(3)

        sigma_ = sigma
        r0_ = r0

        d = norm2(r - r0)

        val = 1._wp/(2._wp*pi*sigma**2)**1.5_wp*exp(-0.5_wp*d**2/sigma**2)
    end subroutine gaussian

    ! Idea solution to the integral \int dV(r') gaussian(r)/4/pi/|r - r'|
    ! the same as the solution to the \nabla^2 u = gaussian
    pure subroutine res_errf(r, val)
        real(wp), intent(in) :: r(3)
        real(wp), intent(out) :: val

        real(wp) :: d
        real(wp) :: sigma_
        real(wp) :: r0_(3)

        sigma_ = sigma
        r0_ = r0

        d = norm2(r - r0)
        val = 0.25_wp/pi/d*erf(d/sqrt(2._wp)/sigma)
    end subroutine res_errf

end program test_integral
