module solver
    use :: typedef, only:wp
    use :: grid_util, only:grid_params
    use :: system_util, only:system
    use :: phys_const, only: &
        cdd => const_cdd, &
        pi => const_pi, &
        one_third => frac_one_third , &
        gLHY => const_gamma_LHY, &
        gcont => const_g_cont

    use, intrinsic :: iso_fortran_env, only: stderr => error_unit

    implicit none
    private

    type, public :: H_GPe
        !!Hamiltionian in the extended Gross-Pitaievskii equation
        real(wp), allocatable :: vdd(:, :, :)
        !!term from dipole dipole interaciton
        real(wp), allocatable :: n(:, :, :)
        !!tabelarized square of the norm of wavefunction
        real(wp), allocatable :: V_ext(:, :, :)
        !!tablearized external trap potential V_ext(x,y,z)

    contains
        procedure, public :: init => init_HGPe
        procedure, public :: apply => apply_HGPe
    end type

    abstract interface
        pure subroutine V_ext0(r, val)
            use :: typedef
            real(wp), intent(in) :: r(3)
            real(wp), intent(out) :: val
        end subroutine V_ext0
    end interface
contains

    !TODO: ensure that the density is in pair with the wavefunction
    !the density refresh is not yet implemented

    subroutine V_dd(H, gr)
        !!compute the dipole diple term of the
        !!extended Gross-Pitajevskii Hamiltonian
        type(H_GPe), intent(inout) :: H
        type(grid_params), intent(in) :: gr

        real(wp), allocatable :: n_copy(:, :, :)
        integer :: i, j, k

        !obnoxious copy, with allocation on assignment
        n_copy = H%n

        !modification in place
        call overrelax_lapl3D(n_copy, gr%dr)

        do k = 2, size(n_copy, 3) - 1
        do j = 1, size(n_copy, 2)
        do i = 1, size(n_copy, 1)
            !central difference for the second derrivative
            !along the dipole axis (z)
            H%vdd(i, j, k) = (n_copy(i, j, k + 1) - 2*n_copy(i, j, k) + n_copy(i, j, k - 1))/gr%dr(3)**2
        end do
        end do
        end do

        !add the term with Kroneker's delta
        !important: here no longer is sys%n equal to n_copy
        H%vdd = H%vdd + one_third*H%n

        !correct scaling constant
        H%vdd = -cdd*H%vdd
    end subroutine V_dd

    subroutine overrelax_Lapl3D(arr, dr)
        !! Perform overrelaxation procedure on for 3D Laplace equation
        real(wp), intent(inout) :: arr(:, :, :)
        !! Grid to be changed in place
        real(wp), intent(in) :: dr(3)
        !! Cell sizes

        integer :: iter
        real(wp) :: vol_inv
        real(wp) :: initial_S
        real(wp) :: final_S

        !TODO: make this smarter. Possibly some convergence-checking criterion
        integer :: N_step = 1000

        !TODO: make this more portable and avoid global state
        real(wp), parameter :: omega = 1.9_wp

        vol_inv = 1._wp/dr(1)/dr(2)/dr(3)

        !set zero boundary conditions
        arr(1, :, :) = 0
        arr(size(arr, 1), :, :) = 0
        arr(:, 1, :) = 0
        arr(:, size(arr, 2), :) = 0
        arr(:, :, 1) = 0
        arr(:, :, size(arr, 3)) = 0

        !initial action functional value
        call calc_action(arr, initial_S)

        !do the relaxation steps
        do iter = 1, N_step
            call step(arr)
        end do

        !determine the final action functional value
        call calc_action(arr, final_S)

        if (final_S > 0.01*initial_S) then
            write (stderr, *) "The relaxation criterion was not met"
            stop 1
        end if
    contains
        pure subroutine step(a)
            real(wp), intent(inout) :: a(:, :, :)
            integer :: i, j, k

            do k = 2, size(a, j) - 1
            do j = 2, size(a, 2) - 1
            do i = 2, size(a, 3) - 1
                a(i, j, k) = (1 - omega)*a(i, j, k) + omega*0.5_wp*vol_inv*( &
                             dr(2)*dr(3)*(a(i - 1, j, k) + a(i + 1, j, k)) + &
                             dr(1)*dr(3)*(a(i, j - 1, k) + a(i, j + 1, k)) + &
                             dr(1)*dr(2)*(a(i, j, k - 1) + a(i, j, k + 1)))
            end do
            end do
            end do
        end subroutine step

        pure subroutine calc_action(a, action)
            real(wp), intent(in) :: a(:, :, :)
            real(wp), intent(out) :: action
            integer :: i, j, k

            action = 0
            do k = 2, size(a, j) - 1
            do j = 2, size(a, 2) - 1
            do i = 2, size(a, 3) - 1
                action = action + abs(0.5_wp*vol_inv*( &
                                      dr(2)*dr(3)*(a(i - 1, j, k) + a(i + 1, j, k)) + &
                                      dr(1)*dr(3)*(a(i, j - 1, k) + a(i, j - 1, k)) + &
                                      dr(1)*dr(2)*(a(i, j, k - 1) + a(i, j, k - 1))) - &
                                      a(i, j, k))
            end do
            end do
            end do
        end subroutine calc_action
    end subroutine overrelax_Lapl3D

    pure subroutine d2_central(arr, dr, diff)
        !!helper subroutine to compute the second derrivative
        !!using central difference
        real(wp), intent(in) :: arr(3)
        !!array containing point of interest and neighbours
        !!along the dimension of intrest
        real(wp), intent(in) :: dr
        !!the size of the grid cell along the dimension
        !!of interest
        real(wp), intent(out) :: diff
        !!computed derrivative value

        diff = (arr(3) - 2*arr(2) + arr(1))/dr**2
    end subroutine d2_central

    subroutine init_HGPe(H, gr, sys, V_ext)
        class(H_GPe), intent(out) :: H
        type(grid_params), intent(in) :: gr
        type(system), intent(in) :: sys
        procedure(V_ext0) :: V_ext

        integer :: i,j,k
        real(wp) :: r(3)

        allocate (H%vdd(gr%dims(1), gr%dims(2), gr%dims(3)))
        allocate (H%V_ext(gr%dims(1), gr%dims(2), gr%dims(3)))
        allocate (H%n(gr%dims(1), gr%dims(2), gr%dims(3)))

        !initialize the external potential
        do i = 1, gr%dims(3)
            r(3) = gr%blc(3) + i*gr%dr(3)
            do j = 1, gr%dims(2)
                r(2) = gr%blc(2) + j*gr%dr(2)
                do k = 1, gr%dims(1)
                    r(1) = gr%blc(1) + k*gr%dr(1)
                    call V_ext(r, H%V_ext(k, j, i))
                end do
            end do
        end do

        !initialize the density
        H%n = sys%psi*conjg(sys%psi)
    end subroutine init_HGPe

    subroutine apply_HGPe(H, grid, sys)
        class(H_GPe), intent(inout) :: H
        type(grid_params), intent(in) :: grid
        type(system), intent(inout) :: sys

        integer :: i, j, k

        !Refresh the stored density
        !additional saety check. Fail if the dimension were to not agree
        !instead of silent reallocation on assgnment
        H%n(:,:,:) = sys%psi*conjg(sys%psi)

        !TODO: find away to possibly avoid obnoxious copy

    end subroutine apply_HGPe
end module solver
