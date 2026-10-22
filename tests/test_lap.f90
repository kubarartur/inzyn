!warning half of this file is vibe coded
program test_lap
    use :: typedef, only:wp
    use :: hamiltonian, only:laplaceian
    use :: grid_util, only:grid_params
    use :: system_util, only:system
    use :: mod_test_lap, only:gaupsi0, psi_analitic
    use, intrinsic :: iso_fortran_env, only: stderr => error_unit
    implicit none

    type(grid_params) :: gr
    type(system) :: sys_num, sys_analitic

    !some random cell size
    real(wp), parameter :: cell_size(3) = [2._wp, 0.7_wp, 3._wp]
    !some random cell numbers
    integer, parameter :: N(3) = [50, 50, 30]
    real(wp), parameter :: tol = 1e-2_wp
    real(wp) :: max_err

    complex(wp), allocatable :: psi_cp(:, :, :)


    type(gaupsi0) :: gausspsi0
    type(psi_analitic) :: psia

    call gr%init(N, cell_size)

    call sys_num%init(gr, gausspsi0, 1)
    allocate (psi_cp, source=sys_num%psi)
    psi_cp = sys_num%psi

    call sys_analitic%init(gr, psia, 1)
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

end program test_lap
