!! module containing definition of computational grid structure
module grid_util
    use :: typedef, only:wp
    use, intrinsic :: iso_fortran_env, only: error_unit

    type, public :: grid
        !!computation grid with parameters:
        real(wp), allocatable :: cell(:, :, :)
        !!actual girdi cells cell(x,y,z)
        real(wp) :: dr(3)
        !!size of the grid cell in each dimension
    contains

        procedure :: init => init_grid
    end type grid

contains

    !> Initialize instance of [[grid]]
    subroutine init_grid(ws, N, dr)
        class(grid), intent(out) :: ws
            !!uninitialized grid
        integer, intent(in) :: N(3)
            !!number of cells in each of the dimensions
        real(wp), intent(in) :: dr(3)
            !!size of the grid cell in each of the dimensions
        allocate (ws%cell(N(1), N(2), N(3)))
        ws%dr = dr
    end subroutine init_grid

end module grid_util
