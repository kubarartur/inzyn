!! module containing definition of computational grid structure
module grid_util
    use :: typedef, only:wp
    use, intrinsic :: iso_fortran_env, only: error_unit

    implicit none
    private

    type, public :: system
        !!system with: grid_step, psi, V_ext
        complex(wp), allocatable :: psi(:, :, :)
        !!tabelarized wavefunction psi(x,y,z)
        real(wp), allocatable :: V_ext(:, :, :)
        !!tablearized external trap potential V_ext(x,y,z)
        real(wp) :: dr(3)
        !!size of the grid cell in each dimension
        real(wp) :: blc(3)
        !!position of the bottom left corner of the grid
    contains
        procedure :: init => init_system
    end type system

    abstract interface
        pure subroutine init_cond_inter_real(r, val)
            use :: typedef
            real(wp), intent(in) :: r(3)
            real(wp), intent(out) :: val
        end subroutine init_cond_inter_real
        pure subroutine init_cond_inter_compl(r, val)
            use :: typedef
            real(wp), intent(in) :: r(3)
            complex(wp), intent(out) :: val
        end subroutine init_cond_inter_compl
    end interface
contains

    !> Initialize instance of [[system]]
    subroutine init_system(ws, N, dr, phi0, V_ext, bottom_left_corner)
        class(system), intent(out) :: ws
            !!uninitialized system
        integer, intent(in) :: N(3)
            !!number of cells in each of the dimensions
        real(wp), intent(in) :: dr(3)
            !!size of the grid cell in each of the dimensions
        procedure(init_cond_inter_compl) :: phi0
            !!initial conditions for the wavefunction
        procedure(init_cond_inter_real) :: V_ext
            !!external trap potential
        real(wp), intent(in), optional :: bottom_left_corner
            !!real position correspoinding to the bottom left corner of the simulation grid

        integer :: i, j, k
        real(wp) :: r(3)

        if(present(bottom_left_corner)) then
            ws%blc = bottom_left_corner
        else
           !if not present make the grid centered
           ws%blc= 0.5_wp*N*dr
        end if

        allocate (ws%psi(N(1), N(2), N(3)))
        allocate (ws%V_ext(N(1), N(2), N(3)))

        ws%dr = dr
        do i = 1, N(3)
            r(3) = k*ws%dr(3)
            do j = 1, N(2)
                r(2) = k*ws%dr(2)
                do k = 1, N(1)
                    r(1) = k*ws%dr(1)
                    call phi0(r, ws%psi(k, j, i))
                    call V_ext(r, ws%V_ext(k, j, i))
                end do
            end do
        end do
    end subroutine init_system

end module grid_util
