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
    real(wp), parameter :: cell_size(3) = [2._wp, 0.7_wp, 3._wp]
    !some random cell numbers
    integer, parameter :: N(3) = [50, 50, 30]
    real(wp), parameter :: tol = 1e-4_wp

    !width and wavevector small enough that O(h^2) FD error stays under tol
    real(wp), parameter :: alpha = 1.0e-4_wp
    real(wp), parameter :: kvec(3) = [2.0e-2_wp, 3.0e-2_wp, 1.5e-2_wp]

    complex(wp), allocatable :: psi_cp(:,:,:)

    call gr%init(N, cell_size)

    call sys_num%init(gr, gauss_psi0)
    allocate(psi_cp, source=sys_num%psi)
    psi_cp = sys_num%psi
    call sys_analitic%init(gr, gauss_laplaceian)

    call laplaceian(gr, sys_num%psi, psi_cp)

    !loop over psi_cp and sys_analitic%psi comparing those two
    block
        integer :: i, j, k
        integer :: imax, jmax, kmax
        real(wp) :: err, max_err

        max_err = 0._wp
        imax = 1
        jmax = 1
        kmax = 1
        do i = 1, gr%dims(3)
        do j = 1, gr%dims(2)
        do k = 1, gr%dims(1)
            err = abs(psi_cp(k, j, i) - sys_analitic%psi(k, j, i))
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

contains

    !initial conditions subroutine conforming to
    !abstract interface ./src/util/system
    !this is reasonably smooth and importantly has a known laplacian
    pure subroutine gauss_psi0(r, val)
        real(wp), intent(in) :: r(3)
        complex(wp), intent(out) :: val

        val = exp(-alpha*sum(r**2))*exp(cmplx(0._wp, sum(kvec*r), kind=wp))
    end subroutine gauss_psi0

    !the analitical laplaceian of those initial conditions
    !abstract interface ./src/util/system
    pure subroutine gauss_laplaceian(r, val)
        real(wp), intent(in) :: r(3)
        complex(wp), intent(out) :: val

        complex(wp) :: psi
        real(wp) :: r2, k2, kr

        r2 = sum(r**2)
        k2 = sum(kvec**2)
        kr = sum(kvec*r)
        call gauss_psi0(r, psi)
        ! nabla^2 [exp(-alpha r^2) exp(i k.r)]
        !   = [4 alpha^2 r^2 - 6 alpha - |k|^2 - 4 i alpha (k.r)] psi
        val = (4._wp*alpha**2*r2 - 6._wp*alpha - k2)*psi
        val = val - cmplx(0._wp, 4._wp*alpha*kr, kind=wp)*psi
    end subroutine gauss_laplaceian
end program test_lap
