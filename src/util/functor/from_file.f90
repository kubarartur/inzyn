module from_file_mod
    use :: typedef, only:wp
    use :: functor, only:compl_functor
    use :: phys_const, only:m_to_a => conv_m_to_a0
    use, intrinsic :: iso_fortran_env, only: error_unit
    implicit none

    type, extends(compl_functor) :: from_file
        real(wp) :: blc(3) = [0._wp, 0._wp, 0._wp]
        real(wp) :: dr(3) = [1._wp, 1._wp, 1._wp]
        complex(wp), allocatable :: values(:, :, :)
    contains
        procedure :: eval => eval_from_file
        procedure :: read_params => read_from_file
    end type

contains
    subroutine eval_from_file(self, r, val)
        !!double gausian initial condition for wavefunction
        !!see implementation for parameters
        class(from_file), intent(in) :: self
        real(wp), intent(in) :: r(3)
        complex(wp), intent(out) :: val

        integer :: indx(3)

        indx = int((r - self%blc)/self%dr)
        val = self%values(indx(1), indx(2), indx(3))

    end subroutine eval_from_file

    subroutine read_from_file(self, fname)
        class(from_file), intent(inout) :: self
        character(len=*), intent(in) :: fname

        integer :: io
        integer :: stat
        integer :: n(3)
        character(len=510) :: line

        integer :: n_part = 20000
        logical :: load_sys = .false.
        character(len=255) :: load_fname = ""
        character(len=255) :: psi0_type = ""
        namelist /system_params/ n_part, load_sys, load_fname, psi0_type

        io = 67
        open (io, file=trim(fname), status="old", action="read")
        read (io, nml=system_params)
        close (io)

        if( .not. load_sys ) then
            write(error_unit, *) "Inside the load from file functor, but the load of system was not requested"
            stop 1
        end if

        open (io, file=trim(load_fname), action="read", status="old", iostat=stat)

        write(error_unit, *) "reading : " , trim(load_fname)
        call read_till("BOTTOM LEFT CORNER")
        read (io, "(A)") line
        read (line, *) self%blc

        call read_till("DIMENSION")
        read (io, "(A)") line
        read (line, *) n

        call read_till("CELL SIZE")
        read (io, "(A)") line
        read (line, *) self%dr

        allocate (self%values(n(1), n(2), n(3)))
        call read_till("ARRAY VALUES")
        read (io, *) self%values

        close (io)
    contains

        subroutine read_till(str)
            character(len=*), intent(in) :: str

            character(len=256) :: line_
            integer :: iostat

            do
                read (io, "(A)", iostat=iostat) line_
                if (iostat /= 0) then
                    close (io)
                    write (error_unit, *) "marker ", str, " not found"
                    stop 1
                end if

                if (index(line_, str) > 0) then
                    exit
                end if
            end do
        end subroutine read_till
    end subroutine read_from_file

end module from_file_mod
