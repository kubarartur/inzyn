!> Basic constants needed for the simulation
!! all constants are prefixed with const_
!! not to clutter the namespace on accident
module phys_const
    use :: typedef, only:wp
    implicit none

    !Unit conversion
    real(wp), parameter :: u_to_me = 1._wp/5.486e-4_wp
    real(wp), parameter :: Hz_to_atm = 4.1356e-12_wp/27211.6_wp
    real(wp), parameter :: energy_renormalization = 1e11
    !real(wp), parameter :: energy_renormalization = 1.e19_wp

    !parameter declaration
    real(wp), parameter :: M = 163.929_wp*u_to_me*energy_renormalization
    real(wp), parameter :: mass = 163.929_wp
    real(wp), parameter :: add = 130 !In a_0 units

    public
    real(wp), parameter :: const_pi = 3.141592653589793_wp

#ifndef NO_OSCILLATOR
    real(wp), parameter :: w(3) = ([60, 60, 120])*Hz_to_atm
    !!external potential trap parameter
#else
    real(wp), parameter :: w(3) = 0._wp
#endif

#ifndef NO_DIPOLE
    real(wp), parameter :: cdd = 12*const_pi*add/M*energy_renormalization*energy_renormalization
#else
    real(wp), parameter :: cdd = 0._wp
#endif

    real(wp), parameter :: epsilon_dd = 1.5
    real(wp), parameter :: a = add/epsilon_dd
contains

end module phys_const
