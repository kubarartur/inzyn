module hamiltonian
    use :: typedef, only:wp
    use :: grid_util, only:grid_params
    use :: system_util, only:system
    use :: phys_const, only: &
        cdd => const_cdd, &
        pi => const_pi, &
        one_third => frac_one_third, &
        gLHY => const_gamma_LHY, &
        gcont => const_g_cont, &
        m => const_m

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

    interface
        module subroutine laplaceian(gr, psi_in, psi_out)
            !!calculate the laplaceian of a complex-valued function psi_in
            !!and write the result to psi_out
            !!function is approximated by values on a grid with gr
            !!for inside points use central approximation of 2nd derrivative
            !!for pointst at the boundary use oneside approximation persevering accuracy
            type(grid_params), intent(in) :: gr
            complex(wp), intent(in) :: psi_in(:, :, :)
            complex(wp), intent(out) :: psi_out(:, :, :)
        end subroutine laplaceian
    end interface
contains

    !TODO: ensure that the density is in pair with the wavefunction
    !the density refresh is not yet implemented

    subroutine apply_HGPe(H, gr, sys)
        class(H_GPe), intent(inout) :: H
        type(grid_params), intent(in) :: gr
        type(system), intent(inout) :: sys

        complex(wp), allocatable :: psi_copy(:, :, :)

        !Refresh the stored density
        !additional saety check. Fail if the dimension were to not agree
        !instead of silent reallocation on assgnment
        H%n(:, :, :) = abs(sys%psi)**2

        !TODO: consider doing something smarter
        !this is expensive maybe though not a bottleneck
        !allocation on assignment
        psi_copy = sys%psi

        call V_dd(H, gr)

        call laplaceian(gr, psi_copy, sys%psi)
        sys%psi = -0.5_wp/m*sys%psi
        sys%psi = sys%psi + (H%V_ext + gcont*H%n + gLHY*H%n*abs(psi_copy) + H%vdd)*psi_copy
    end subroutine apply_HGPe

    subroutine V_dd(H, gr)
        !!compute the dipole diple term of the
        !!extended Gross-Pitajevskii Hamiltonian
        type(H_GPe), intent(inout) :: H
        type(grid_params), intent(in) :: gr

        real(wp), allocatable :: n_copy(:, :, :)
        integer :: i, j, k

        !obnoxious copy, with allocation on assignment
        allocate (n_copy, source=H%n)
        n_copy = H%n

        !modification in place
        call overrelax_lapl3D(n_copy, H%n, gr%dr)

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
        !important: here no longer is H%n equal to n_copy
        H%vdd = H%vdd + one_third*H%n

        !correct scaling constant
        H%vdd = -cdd*H%vdd
    end subroutine V_dd

    subroutine overrelax_Lapl3D(arr, rhs, dr)
        !! Perform overrelaxation procedure for 3D Laplace equation
        real(wp), intent(inout) :: arr(:, :, :)
        !! Grid to be changed in place
        real(wp), intent(in) :: rhs(:, :, :)
        !! The right hand side of the laplace's equation to be solved
        real(wp), intent(in) :: dr(3)
        !! Cell sizes

        integer :: iter
        real(wp) :: initial_S
        real(wp) :: final_S

        !TODO: figure out a smarter name for this variable
        real(wp) :: multip

        !TODO: make this smarter. Possibly some convergence-checking criterion
        integer :: N_step = 1000

        !TODO: make this more portable
        real(wp), parameter :: omega = 1.9_wp

        multip = 0.5_wp/(1/dr(1)**2 + 1/dr(2)**2 + 1/dr(3)**2)
        !set zero boundary conditions
        arr(1, :, :) = 0
        arr(size(arr, 1), :, :) = 0
        arr(:, 1, :) = 0
        arr(:, size(arr, 2), :) = 0
        arr(:, :, 1) = 0
        arr(:, :, size(arr, 3)) = 0

        !initial action functional value
        !TODO: Verify that calc action indeed works correctly
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
        subroutine step(a)
            real(wp), intent(inout) :: a(:, :, :)
            integer :: i, j, k

            do k = 2, size(a, 3) - 1
            do j = 2, size(a, 2) - 1
            do i = 2, size(a, 1) - 1
                a(i, j, k) = (1 - omega)*a(i, j, k) + omega*multip*( &
                             (a(i - 1, j, k) + a(i + 1, j, k))/dr(1)**2 + &
                             (a(i, j - 1, k) + a(i, j + 1, k))/dr(2)**2 + &
                             (a(i, j, k - 1) + a(i, j, k + 1))/dr(3)**2 - &
                             rhs(i, j, k))
            end do
            end do
            end do
        end subroutine step

        pure subroutine calc_action(a, action)
            real(wp), intent(in) :: a(:, :, :)
            !!current approximation
            real(wp), intent(out) :: action
            !!action functional value
            integer :: i, j, k

            real(wp) :: mult(3)
            real(wp) :: vol
            vol = dr(1)*dr(2)*dr(3)
            mult(1) = dr(2)*dr(3)/dr(1)
            mult(2) = dr(1)*dr(3)/dr(2)
            mult(3) = dr(1)*dr(2)/dr(3)

            action = 0
            do k = 2, size(a, 3) - 1
            do j = 2, size(a, 2) - 1
            do i = 2, size(a, 1) - 1
                action = action + 0.125_wp*( &
                         mult(1)*(a(i + 1, j, k) - a(i - 1, j, k))**2 + &
                         mult(2)*(a(i, j + 1, k) - a(i, j - 1, k))**2 + &
                         mult(3)*(a(i, j, k + 1) - a(i, j, k - 1))**2 &
                         ) + rhs(i, j, k)*a(i, j, k)*vol
            end do
            end do
            end do
        end subroutine calc_action
    end subroutine overrelax_Lapl3D

    subroutine init_HGPe(H, gr, sys, V_ext)
        class(H_GPe), intent(out) :: H
        type(grid_params), intent(in) :: gr
        type(system), intent(in) :: sys
        procedure(V_ext0) :: V_ext

        integer :: i, j, k
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
        H%n = abs(sys%psi)**2
    end subroutine init_HGPe

end module hamiltonian
