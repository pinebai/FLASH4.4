!!****if* source/physics/Gravity/GravityMain/Poisson/BHTree/grv_bhAccShort
!!
!! NAME
!!
!!  grv_bhAccShort
!!
!!
!! SYNOPSIS
!!
!!   real(0:12) ef = grv_bhAccShort(
!!                       real(in) :: ewald_alpha,
!!                       real(in) :: xni,
!!                       real(in) :: yni,
!!                       real(in) :: zni
!!                   )
!!
!! DESCRIPTION
!!
!!   Calculates potential, acceleration and partial derivatives of acceleration 
!!   for the Ewald field for short-range contributions.
!!
!! ARGUMENTS
!!
!!   ewald_alpha  : constant controlling convergence rate
!!   xni          : relative distance between target point and point mass in x direction
!!   yni          : relative distance between target point and point mass in y direction
!!   zni          : relative distance between target point and point mass in z direction
!!
!! RESULT
!!
!!   Returns an array with 13 real numbers. 
!!   Index 0 is a value of potential
!!   indices 1 to 3 are values of acceleration
!!   indices 4 to 12 are values of partial derivatives of acceleration
!!   Note that three derivatives (at 7, 10 and 11) are not evaluated because 
!!   of the symmetry.
!!   
!!
!!***

#include "Flash.h"

function grv_bhAccShort(ewald_alpha,xni,yni,zni)

#if defined FLASH_USE_SPECFUN
  use grv_bhInterface, ONLY : grv_derfc
#endif
  implicit none
#include "constants.h"

  real, intent(IN) :: ewald_alpha,xni,yni,zni
  real :: grv_bhAccShort(0:12)

  real, parameter :: pi = PI

  real :: rni, rni2, erfcar
  real :: a1, a2, a3, a4

#ifdef USER_ERFC
#define ERFC USER_ERFC
  real, external :: USER_ERFC
#elif defined(FLASH_USE_SPECFUN)
#define ERFC grv_derfc
#else
#define ERFC erfc
  intrinsic erfc
#endif

  rni2 = xni*xni + yni*yni + zni*zni
  rni = sqrt(rni2)
  erfcar = ERFC(ewald_alpha*rni)
  a1 = 2.0D0*ewald_alpha*exp(-rni2*ewald_alpha**2)/(sqrt(pi)*rni2)
  a2 = erfcar/rni**3
  a3 = erfcar/rni**5
  a4 = 3.0D0/rni2 + 2.0D0*ewald_alpha*ewald_alpha

  ! potential
  grv_bhAccShort(0) = erfcar / rni
  ! acceleration
  grv_bhAccShort(1) = (a1 + a2)*xni
  grv_bhAccShort(2) = (a1 + a2)*yni
  grv_bhAccShort(3) = (a1 + a2)*zni
  ! partial derivatives of acceleration
  grv_bhAccShort(4) = a2 - 3*(xni**2)*a3 + a1*(1.0D0 - a4*xni**2)
  grv_bhAccShort(5) = -xni*yni*(3*a3 + a1*a4)
  grv_bhAccShort(6) = -xni*zni*(3*a3 + a1*a4)
  grv_bhAccShort(8) = a2 - 3*(yni**2)*a3 + a1*(1.0D0 - a4*yni**2)
  grv_bhAccShort(9) = -yni*zni*(3*a3 + a1*a4)
  grv_bhAccShort(12) = a2 - 3*(zni**2)*a3 + a1*(1.0D0 - a4*zni**2)

end


