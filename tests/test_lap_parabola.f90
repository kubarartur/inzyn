!warning half of this file is vibe coded
program test_lap
    use :: typedef, only:wp
    use :: hamiltonian, only:laplaceian
    use :: grid_util, only:grid_params
    use :: system_util, only:system
    use, intrinsic :: iso_fortran_env, only: stderr => error_unit
    implicit none

    type(grid_params) :: gr
    type(system) :: sys_num, sys_analitic

    !some random cell size
    real(wp), parameter :: cell_size(3) = [90._wp, 70._wp, 300._wp]
    !some random cell numbers
    integer, parameter :: N(3) = [50, 50, 30]
    real(wp), parameter :: tol = 1e-2_wp
    real(wp) :: max_err

    complex(wp), allocatable :: psi_cp(:, :, :)

    integer :: ii
    real(wp), target :: coefs_flat(18)
    real(wp), pointer :: coefs(:, :, :)
    !! r/i, dim abc
    integer :: seed_size
    integer, allocatable :: my_seed(:)

    call random_seed(size=seed_size)
    allocate (my_seed(seed_size))
    my_seed = 1
    call random_seed(put=my_seed)

    do ii = 1, 18
        call random_number(coefs_flat(ii))
        coefs_flat = (coefs_flat - 0.5_wp)*21.37_wp
    end do
    coefs(1:2, 1:3, 1:3) => coefs_flat

    call gr%init(N, cell_size)

    call sys_num%init(gr, parabola)
    allocate (psi_cp, source=sys_num%psi)
    psi_cp = sys_num%psi

    call sys_analitic%init(gr, parabola_diff)
    call laplaceian(gr, sys_num%psi, psi_cp)

    !loop over psi_cp and sys_analitic%psi comparing those two
    block
        integer :: i, j, k
        integer :: imax, jmax, kmax
        real(wp) :: err

        max_err = 0._wp
        imax = 1
        jmax = 1
        kmax = 1
        do i = 1, gr%dims(3)
        do j = 1, gr%dims(2)
        do k = 1, gr%dims(1)
            err = abs(psi_cp(k, j, i) - sys_analitic%psi(k, j, i))/ &
                  abs(sys_analitic%psi(k, j, i) + 1e-8_wp)
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
            write (stderr, *) "laplacian mismatch at index: ", kmax, jmax, imax
            write (stderr, *) "numerical:  ", psi_cp(kmax, jmax, imax)
            write (stderr, *) "analytical: ", sys_analitic%psi(kmax, jmax, imax)
            write (stderr, *) "The error of: ", max_err
            write (stderr, *) "exceeded the tolerance of: ", tol
            stop 1
        end if
    end block
    write (*, *) "Passed with max error of: ", max_err, "for tol of: ", tol

contains

    pure subroutine parabola(r, val)
        real(wp), intent(in) :: r(3)
        complex(wp), intent(out) :: val

        !TODO: completely inefficient index order
        val = complex(sum(coefs(1, 1, :)*[r(1)**2, r(1), 1._wp]), sum(coefs(2, 1, :)*[r(1)**2, r(1), 1._wp]))
        val = val + complex(sum(coefs(1, 2, :)*[r(2)**2, r(2), 1._wp]), sum(coefs(2, 2, :)*[r(2)**2, r(2), 1._wp]))
        val = val + complex(sum(coefs(1, 3, :)*[r(3)**2, r(3), 1._wp]), sum(coefs(2, 3, :)*[r(3)**2, r(3), 1._wp]))
    end subroutine parabola

    !should be equal to the laplacian of parabola
    pure subroutine parabola_diff(r, val)
        real(wp), intent(in) :: r(3)
        complex(wp), intent(out) :: val

        !too keep the compiler silent about unused
        if (.false.) then
            val = complex(sum(r), 0._wp)
        end if

        val = complex(sum(coefs(1, 1, :)*[2._wp, 0._wp, 0._wp]), sum(coefs(2, 1, :)*[2._wp, 0._wp, 0._wp]))
        val = val + complex(sum(coefs(1, 2, :)*[2._wp, 0._wp, 0._wp]), sum(coefs(2, 2, :)*[2._wp, 0._wp, 0._wp]))
        val = val + complex(sum(coefs(1, 3, :)*[2._wp, 0._wp, 0._wp]), sum(coefs(2, 3, :)*[2._wp, 0._wp, 0._wp]))
    end subroutine parabola_diff
end program test_lap
