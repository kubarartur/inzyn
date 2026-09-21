!> Basic constants needed for the simulation
module phys_const
   use my_kind
   implicit none
   !Unit conversion
   real(wp), parameter :: u_to_me = 1._wp/5.486e-4_wp
   real(wp), parameter :: Hz_to_atm = 4.1356e-12_wp/27211.6_wp
   real(wp), parameter :: pi = 3.141592653589793_wp
   real(wp), parameter :: energy_renormalization = 1e11
   !real(wp), parameter :: energy_renormalization = 1.e19_wp

   !parameter declaration
   real(wp), parameter :: M = 163.929_wp*u_to_me*energy_renormalization

#ifndef NO_OSCILLATOR
   real(wp), parameter :: w(3) = ([60, 60, 120])*Hz_to_atm
#else
   real(wp), parameter :: w(3) = 0._wp !([60 , 60, 120])*Hz_to_atm
#endif

   real(wp), parameter :: mass = 163.929_wp
   real(wp), parameter :: add = 130 !In a_0 units
   ! some wierdly indented comment
#ifndef NO_DIPOLE
   real(wp), parameter :: cdd = 12*pi*add/M*energy_renormalization*energy_renormalization
#else
   real(wp), parameter :: cdd = 0._wp
#endif

   real(wp), parameter :: epsilon_dd = 1.5
   real(wp), parameter :: a = add/epsilon_dd

   !Repulsive potential
#ifndef NO_CONTACT
   real(wp), parameter :: len_eps = 8e-11_wp*energy_renormalization
    !! $\sigma$ parameter of the lenard jones potential term $\frcac{1}{d^6}$
#else
   real(wp), parameter :: len_eps = 0._wp
    !! $\sigma$ parameter of the lenard jones potential term $\frcac{1}{d^6}$
#endif

   !fractions
   real(wp), parameter :: one_third = 0.333333333333333333333333_wp
   real(wp), parameter :: one_sixth = 0.166666666666666666666666_wp
contains

end module phys_const
