!> Basic constants needed for the simulation
!! all constants are prefixed with const_
!! not to clutter the namespace on accident
module phys_const
    use :: typedef, only:wp
    implicit none
    private

    !Unit conversion
    real(wp), parameter :: u_to_me = 1._wp/5.486e-4_wp
    real(wp), parameter :: Hz_to_atm = 4.1356e-12_wp/27211.6_wp
    real(wp), parameter, public :: conv_m_to_a0= 5.29177210544e-11_wp

    !parameter declaration
    real(wp), parameter :: add = 130 !In a_0 units

    real(wp), parameter, public :: const_m = 163.929_wp*u_to_me

    real(wp), parameter, public :: const_pi = 3.141592653589793_wp

#ifndef NO_OSCILLATOR
    real(wp), parameter, public  :: const_w(3) = ([60, 60, 120])*Hz_to_atm
    !!external potential trap parameter
#else
    real(wp), parameter, public  :: const_w(3) = 0._wp
#endif

#ifndef NO_DIPOLE
    real(wp), parameter, public  :: const_cdd = 12*const_pi*add/const_m
#else
    real(wp), parameter, public  :: cosnt_cdd = 0._wp
#endif

    real(wp), parameter :: epsilon_dd = 1.5
    real(wp), parameter :: a = add/epsilon_dd
contains

end module phys_const
