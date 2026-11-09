module from_file_mod
    use :: typedef, only:wp
    use :: functor, only:compl_functor
    use :: phys_const, only:m_to_a => conv_m_to_a0
    use :: arr_iout, only:load_arr
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

        indx = int((r - self%blc)/self%dr + 0.5_wp)

        !TODO: add index out of bounds checking
        val = self%values(indx(1), indx(2), indx(3))
    end subroutine eval_from_file

    subroutine read_from_file(self, fname)
        class(from_file), intent(inout) :: self
        character(len=*), intent(in) :: fname

        integer :: io

        integer :: n_part = 20000
        logical :: load_sys = .false.
        character(len=255) :: load_fname = ""
        character(len=255) :: psi0_type = ""

        real(wp) :: bottom_left_corner(3) = [0._wp, 0._wp, 0._wp]
        real(wp) :: cell_size(3) = [0._wp, 0._wp, 0._wp]

        namelist /system_params/ n_part, load_sys, load_fname, psi0_type
        namelist /spacial_params/ bottom_left_corner, cell_size

        io = 67
        open (io, file=trim(fname), status="old", action="read")
        read (io, nml=system_params)
        close (io)

        if (.not. load_sys) then
            write (error_unit, *) "Inside the load from file functor, but the load of system was not requested"
            stop 1
        end if

        open (io, file=trim(load_fname)//"_spatial.nml", status="old", action="read")
        read (io, nml=spacial_params)
        close (io)

        self%blc = bottom_left_corner
        self%dr = cell_size

        write(*,*) "here"
        self%blc = bottom_left_corner
        self%dr = cell_size

        call load_arr(trim(load_fname), self%values)
    end subroutine read_from_file

end module from_file_mod
