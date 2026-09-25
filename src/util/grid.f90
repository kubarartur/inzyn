!! module containing definition of computational grid structure
module grid_util
    use :: typedef, only:wp
    use, intrinsic :: iso_fortran_env, only: error_unit

    implicit none
    private

    type, public :: grid_params
        real(wp) :: dr(3)
        !!size of the grid cell in each dimension
        real(wp) :: blc(3)
        !!position of the bottom left corner of the grid
        real(wp), allocatable :: dist(:, :, :)
        !!table of distances based on the relative offset
        !!to persevere logic lower bound is 0
        integer :: dims(3)
    contains

        procedure :: init => init_grid
        !!initialize the grid

    end type grid_params

contains

    !> Initialize instance of [[system]]
    subroutine init_grid(gr, N, dr, bottom_left_corner)
        class(grid_params), intent(out) :: gr
            !!uninitialized grid_params instance
        integer, intent(in) :: N(3)
            !!number of cells in each of the dimensions
        real(wp), intent(in) :: dr(3)
            !!size of the grid cell in each dimension
        real(wp), intent(in), optional :: bottom_left_corner
            !!real position correspoinding to the bottom left corner of the simulation grid
            !!to make the arithme tic easy this is the center of imaginary cell(0,0,0)
            !! r{cell(i,j,k)}= blc + [i,j,k]*dr

        integer :: i, j, k
        real(wp) :: r(3)

        gr%dr = dr
        gr%dims = N
        if (present(bottom_left_corner)) then
            gr%blc = bottom_left_corner
        else
            !if not present make the grid centered
            gr%blc = -0.5_wp*(N + 1)*dr
        end if

        allocate (gr%dist(0:N(1) - 1, &
                          0:N(2) - 1, &
                          0:N(3) - 1))

        !tabelarize the distances
        do i = 0, N(3) - 1
            r(3) = i*gr%dr(3)
            do j = 0, N(2) - 1
                r(2) = j*gr%dr(2)
                do k = 0, N(1) - 1
                    r(1) = k*gr%dr(1)
                    gr%dist(k, j, i) = norm2(r)
                end do
            end do
        end do

    end subroutine init_grid

end module grid_util
