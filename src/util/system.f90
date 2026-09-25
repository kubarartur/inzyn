module system_util
    use :: typedef, only:wp
    use :: grid_util, only:grid_params
    implicit none
    private

    type, public :: system
        !!system with: psi
        complex(wp), allocatable :: psi(:, :, :)
    contains
        procedure :: init => init_system
    end type system

    abstract interface
        pure subroutine init_cond_psi0(r, val)
            use :: typedef
            real(wp), intent(in) :: r(3)
            complex(wp), intent(out) :: val
        end subroutine init_cond_psi0
    end interface

contains
    subroutine init_system(sys, gr, psi0)
        class(system), intent(out) :: sys
        type(grid_params), intent(in) :: gr
        procedure(init_cond_psi0) :: psi0

        integer :: i, j, k
        real(wp) :: r(3)

        allocate (sys%psi(gr%dims(1), gr%dims(2), gr%dims(3)))

        !initialize the wavefunction
        do i = 1, gr%dims(3)
            r(3) = gr%blc(3) + i*gr%dr(3)
            do j = 1, gr%dims(2)
                r(2) = gr%blc(2) + j*gr%dr(2)
                do k = 1, gr%dims(1)
                    r(1) = gr%blc(1) + k*gr%dr(1)
                    call psi0(r, sys%psi(k, j, i))
                end do
            end do
        end do
    end subroutine init_system

end module
