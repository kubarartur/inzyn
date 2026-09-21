!> wp precission definition
!! infinity and signaling nan utility functions
module typedef
    use, intrinsic :: iso_fortran_env, only: wp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_signaling_nan, ieee_positive_inf
    implicit none
contains

    pure function get_snan() result(nan)
      !! get signaling nan of real working precison type
        real(wp) :: nan
          !! working precision signaling nan scalar
        nan = ieee_value(1.0_wp, ieee_signaling_nan)
    end function get_snan

    pure function get_inf() result(inf)
     !! get infiinity of real working precion type
        real(wp) :: inf
          !! working precision infinity scalar
        inf = ieee_value(1.0_wp, ieee_positive_inf)
    end function get_inf

end module typedef
