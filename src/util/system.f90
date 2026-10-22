module system_util
    use :: typedef, only:wp
    use :: grid_util, only:grid_params
    use :: initial_conditions, only:compl_functor
    implicit none
    private

    type, public :: system
        !!system with: psi
        complex(wp), allocatable :: psi(:, :, :)
        integer :: ncz
        !!integer number of particles
        !TODO: consider weather to keep this or not
        real(wp) :: nczr
        !! real copy of the number of particles to compute sqrt
    contains
        procedure :: init => init_system
        procedure :: copy => copy_system
    end type system


contains
    subroutine init_system(sys, gr, psi0, ncz)
        class(system), intent(out) :: sys
        type(grid_params), intent(in) :: gr
        class(compl_functor) :: psi0
        integer, intent(in) :: ncz

        integer :: i, j, k
        real(wp) :: r(3)

        sys%ncz = ncz
        sys%nczr = real(sys%ncz)
        allocate (sys%psi(gr%dims(1), gr%dims(2), gr%dims(3)))

        !initialize the wavefunction
        do i = 1, gr%dims(3)
            r(3) = gr%blc(3) + i*gr%dr(3)
            do j = 1, gr%dims(2)
                r(2) = gr%blc(2) + j*gr%dr(2)
                do k = 1, gr%dims(1)
                    r(1) = gr%blc(1) + k*gr%dr(1)
                    call psi0%eval(r, sys%psi(k, j, i))
                end do
            end do
        end do
    end subroutine init_system

    subroutine copy_system(sys, sys2)
        class(system), intent(out) :: sys
        type(system), intent(in) :: sys2
        !TODO: possifbly fragile as one argument is polymorphic and other not
        ! at the time of writing had no better idea

       allocate(sys%psi, source=sys2%psi)
        sys%psi(:,:,:) = sys2%psi
    end subroutine copy_system
end module
