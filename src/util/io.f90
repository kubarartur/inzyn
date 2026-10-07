module io
    use :: typedef, only:wp
    use :: grid_util, only:grid_params
    implicit none
    private

    public :: write_arr
    interface write_arr
        module procedure write_arr_real3d
        module procedure write_arr_compl3d
        module procedure write_arr_compl1d
    end interface
contains

    subroutine write_arr_real3d(arr, gr, descr, fname)
        real(wp), intent(in) :: arr(:, :, :)
        type(grid_params), intent(in) :: gr
        character(len=*), intent(in) :: descr
        !!spaces here will break my parser :(
        character(len=*), intent(in) :: fname

        integer :: io
        logical :: ex
        io = 123

        inquire (file=fname, exist=ex)
        if (ex) then
            open (newunit=io, file=fname, status="replace", action="write")
        else
            open (newunit=io, file=fname, status="new", action="write")
        end if

        write (io, *) "DESCRIPTION: ", descr
        write (io, *) "BOTTOM LEFT CORNER OF THE SIMULATIO GRID"
        write (io, *) gr%blc
        write (io, *) "DIMENSION OF THE SIMULATION GRID"
        write (io, *) gr%dims
        write (io, *) "CELL SIZE"
        write (io, *) gr%dr
        write (io, *) "ARRAY VALUES (REAL): "
        write (io, *) arr

        close (io)
    end subroutine write_arr_real3d

    subroutine write_arr_compl3d(arr, gr, descr, fname)
        complex(wp), intent(in) :: arr(:, :, :)
        type(grid_params), intent(in) :: gr
        character(len=*), intent(in) :: descr
        !!spaces here will break my parser :(
        character(len=*), intent(in) :: fname

        integer :: io
        logical :: ex
        io = 123

        inquire (file=fname, exist=ex)
        if (ex) then
            open (newunit=io, file=fname, status="replace", action="write")
        else
            open (newunit=io, file=fname, status="new", action="write")
        end if

        write (io, *) "DESCRIPTION: ", descr
        write (io, *) "BOTTOM LEFT CORNER OF THE SIMULATIO GRID"
        write (io, *) gr%blc
        write (io, *) "DIMENSION OF THE SIMULATION GRID"
        write (io, *) gr%dims
        write (io, *) "CELL SIZE"
        write (io, *) gr%dr
        write (io, *) "ARRAY VALUES (COMPLEX): "
        write (io, *) arr

        close (io)
    end subroutine write_arr_compl3d

    subroutine write_arr_compl1d(arr, extent, descr, fname)
        complex(wp), intent(in) :: arr(:)
        real(wp), intent(in) :: extent(2)
        character(len=*), intent(in) :: descr
        !!spaces here will break my parser :(
        character(len=*), intent(in) :: fname

        integer :: io
        logical :: ex
        io = 123

        inquire (file=fname, exist=ex)
        if (ex) then
            open (newunit=io, file=fname, status="replace", action="write")
        else
            open (newunit=io, file=fname, status="new", action="write")
        end if

        write (io, *) "DESCRIPTION: ", descr
        write (io, *) "EXTENT: "
        write (io, *) extent
        write (io, *) "DIMENSION"
        write (io, *) shape(arr)
        write (io, *) "STEP SIZE"
        write (io, *) (extent(2) - extent(1))/size(arr, 1)
        write (io, *) "ARRAY VALUES (COMPLEX): "
        write (io, *) arr

        close (io)
    end subroutine write_arr_compl1d
end module io
