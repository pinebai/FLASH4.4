!!****if* source/physics/Gravity/GravityMain/Poisson/BHTree/grv_intCoef1P2I
!!
!! NAME
!!
!!  grv_intCoef1P2I
!!
!!
!! SYNOPSIS
!!
!!    real = grv_intCoef1P2I(
!!                           real(in)    :: x,
!!                           integer(in) :: l,
!!                           real(in)    :: eta,
!!                           real(in)    :: dzeta
!!        )
!!
!!    real = grv_intCoef1P2I_der(
!!                           real(in)    :: x,
!!                           integer(in) :: l,
!!                           real(in)    :: eta,
!!                           real(in)    :: dzeta
!!        )
!!
!!    real = grv_intCoef1P2I_der2(
!!                           real(in)    :: x,
!!                           integer(in) :: l,
!!                           real(in)    :: eta,
!!                           real(in)    :: dzeta
!!        )
!!
!!
!! DESCRIPTION
!!
!!   grv_intCoef1P2I evaluates functions which are numerically
!!   integrated by routine grv_intSimpson. They are relevant to long-range 
!!   interactions in the case with periodic boundaries in 1 direction
!!   and isolated boundaries in the other two directions.
!!   
!!
!! ARGUMENTS
!!
!!  l            - wavenumber in the periodic direction
!!  eta          - normalised distance from the cylinder axis
!!  dzeta        - parameter related to the Ewald method. 
!!                 It is defined in Wunsch et al. 2015 (section 2.2.3)
!!  x            - normalised wavevector in the isolated directions
!!
!! RESULT
!!
!!  Integrated functions in the expression for long-range interaction 
!!  in the case with mixed boundary conditions with one periodic direction.
!!
!! NOTES
!!
!!***

#include "Flash.h"

! integrated functions
real function grv_coef1P2I(x,l,eta,dzeta)
#if defined FLASH_USE_SPECFUN
  use grv_bhInterface, ONLY : grv_besj0
#elif defined(__INTEL_COMPILER)
  use IFPORT
#endif

  implicit none

  real,intent(in) :: eta,dzeta,x
  integer,intent(in) :: l

#if defined USER_BESSEL_J0
#define BESSEL_J0 USER_BESSEL_J0
  real, external :: USER_BESSEL_J0
#elif defined(FLASH_USE_SPECFUN)
#define BESSEL_J0 grv_besj0
#else
#define BESSEL_J0 bessel_j0
  intrinsic bessel_j0
#endif

  if (l.eq.0) then
    grv_coef1P2I=(BESSEL_J0(x*eta)-1.0)*exp(-dzeta*x*x)/x
  else
    grv_coef1P2I=x*BESSEL_J0(x*eta)*exp(-dzeta*x*x)/(x**2+l**2)
  endif

end

real function grv_coef1P2I_der(x,l,eta,dzeta)
#if defined FLASH_USE_SPECFUN
  use grv_bhInterface, ONLY : grv_besj1
#elif defined(__INTEL_COMPILER)
  use IFPORT
#endif

  implicit none

  real,intent(in) :: eta,dzeta,x
  integer,intent(in) :: l

#if defined USER_BESSEL_J1
#define BESSEL_J1 USER_BESSEL_J1
  real, external :: USER_BESSEL_J1
#elif defined(FLASH_USE_SPECFUN)
#define BESSEL_J1 grv_besj1
#else
#define BESSEL_J1 bessel_j1
  intrinsic bessel_j1
#endif

  grv_coef1P2I_der=BESSEL_J1(x*eta)*exp(-dzeta*x*x)*(x**2)/(x**2+l**2)

end

real function grv_coef1P2I_der2(x,l,eta,dzeta)

#if defined FLASH_USE_SPECFUN
  use grv_bhInterface, ONLY : grv_besj0, grv_besj1
#elif defined(__INTEL_COMPILER)
  use IFPORT
#endif

  implicit none

#if defined USER_BESSEL_J0
#define BESSEL_J0 USER_BESSEL_J0
  real, external :: USER_BESSEL_J0
#elif defined(FLASH_USE_SPECFUN)
#define BESSEL_J0 grv_besj0
#else
#define BESSEL_J0 bessel_j0
  intrinsic bessel_j0
#endif

#if defined USER_BESSEL_J1
#define BESSEL_J1 USER_BESSEL_J1
  real, external :: USER_BESSEL_J1
#elif defined(FLASH_USE_SPECFUN)
#define BESSEL_J1 grv_besj1
#else
#define BESSEL_J1 bessel_j1
 intrinsic bessel_j1
#endif

 real,intent(in) :: eta,dzeta,x
 integer,intent(in) :: l


  grv_coef1P2I_der2 = (BESSEL_J0(x*eta) - BESSEL_J1(x*eta)/(x*eta))*exp(-dzeta*x*x)* &
  &                   (x**3)/(x**2+l**2)

end



