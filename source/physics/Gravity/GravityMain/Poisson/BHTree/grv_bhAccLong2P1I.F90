!!****if* source/physics/Gravity/GravityMain/Poisson/BHTree/grv_bhAccLong2P1I
!!
!! NAME
!!
!!  grv_bhAccLong2P1I
!!
!!
!! SYNOPSIS
!!
!!   real(0:12) ef = grv_bhAccLong2P1I(
!!                       real(in) :: Linv,
!!                       real(in) :: ewald_dzeta,
!!                       real(in) :: ratio_pinv,
!!                       integer(in) :: hi,
!!                       integer(in) :: hj,
!!                       real(in) :: x,
!!                       real(in) :: y,
!!                       real(in) :: z
!!                   )
!!
!! DESCRIPTION
!!
!!   Calculates potential, acceleration and partial derivatives of acceleration 
!!   for the Ewald field for long-range contributions in the case of mixed boundary 
!!   conditions - periodic in 2 directions and isolated in the third direction.
!!
!! ARGUMENTS
!!
!!   Linv         : inverted size of the computational domain in one of two periodic
!!                  directions
!!   ewald_dzeta  : coefficient dzeta needed to calculate the Ewald field
!!   ratio_pinv   : ratio between sizes of computational domain in the periodic directions
!!   hi           : wavenumber in one of the periodic directions
!!   hj           : wavenumber in the other periodic direction
!!   x            : x-coordinate of the point where the Ewald field is determined
!!   y            : y-coordinate of the point where the Ewald field is determined
!!   z            : z-coordinate of the point where the Ewald field is determined
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

function grv_bhAccLong2P1I(Linv, ewald_dzeta, ratio_pinv, hi, hj, x, y, z)

#if defined FLASH_USE_SPECFUN
  use grv_bhInterface, ONLY : grv_derfc, grv_derfcx, grv_derf
#elif defined(__INTEL_COMPILER)
  use IFPORT
#endif
  implicit none
#include "constants.h"
  real, parameter :: pi = PI

  real, intent(IN) :: Linv, ewald_dzeta, ratio_pinv, x, y, z
  integer, intent(IN) :: hi, hj
  real :: grv_bhAccLong2P1I(0:12)

  real :: hjr, ewc1, ewc2, ewc3, ewc4, ewc5, ewc6, ewc7
  real :: ewald_betas, ewald_betac, ewald_gamma, Linv2
  real :: grav_integral, der_grav_integral, der2_grav_integral

#ifdef USER_ERFC
#define ERFC USER_ERFC
  real, external :: USER_ERFC
#elif defined(FLASH_USE_SPECFUN)
#define ERFC grv_derfc
#else
#define ERFC erfc
  intrinsic erfc
#endif

#ifdef USER_ERFC_SCALED
#define ERFC_SCALED USER_ERFC_SCALED
  real, external :: USER_ERFC_SCALED
#elif defined(FLASH_USE_SPECFUN)
#define ERFC_SCALED grv_derfcx
#else
#define ERFC_SCALED erfc_scaled
  intrinsic erfc_scaled
#endif

#ifdef USER_ERF
#define ERF USER_ERF
  real, external :: USER_ERF
#elif defined(FLASH_USE_SPECFUN)
#define ERF grv_derf
#else
#define ERF erf
  intrinsic erf
#endif

! If compilation fails in one of the preceding lines, or there are other problems
! related to ERFC, see file README.erfc in the BHTree/Wunsch implementation
! directory!

  Linv2 = 2.0*Linv*ratio_pinv/(pi)
  ewald_gamma = 2*pi*z*Linv

!  constants which simplify equations below
  ewc1 = 0.25*pi
  ewc2 = 1.0/sqrt(ewald_dzeta)
  ewc3 = 0.5*ewc2
  ewc4 = 2.0*sqrt(ewald_dzeta/pi)
  ewc5 = 0.5*sqrt(pi*ewald_dzeta)
  ewc6 = 0.25/ewald_dzeta

  if ((hi == 0).and.(hj == 0)) then
    grav_integral = 2.0*(ewc1*(ewald_gamma*ERF(ewc3*ewald_gamma)+ &
    & ewc4*exp(-ewc6*ewald_gamma*ewald_gamma)))
    der_grav_integral = pi*ERF(ewald_gamma*ewc3)
    der2_grav_integral = 2*Linv*(pi**2)*exp(-ewc6*ewald_gamma*ewald_gamma)/(sqrt(pi*ewald_dzeta))

! only three terms are non-zero
    grv_bhAccLong2P1I(:) = 0.0D0
    grv_bhAccLong2P1I(0) = -grav_integral*Linv2
    grv_bhAccLong2P1I(3) =  2*(Linv**2)*ratio_pinv*der_grav_integral
    grv_bhAccLong2P1I(12) =  2*(Linv**2)*ratio_pinv*der2_grav_integral
  else
    hjr = ratio_pinv*real(hj)
    ewald_betas = sin((2*pi)*(hi*x+hjr*y)*Linv)
    ewald_betac = cos((2*pi)*(hi*x+hjr*y)*Linv)

    ewc7 = sqrt(real(hi**2 + hjr**2))
    ! analytical evaluation of appropriate integrals
    ! (in order to avoid behaviour like infty*zero we express
    ! the function with ERFC_SCALED function where needed)

    grav_integral = ewc1*(ERFC_SCALED(ewc2*(ewc7*ewald_dzeta+0.5*ewald_gamma)) &
    &               * exp(-ewald_gamma**2/(4.0D0*ewald_dzeta)-ewald_dzeta*ewc7*ewc7) &
    &               + ERFC(ewc2*(ewc7*ewald_dzeta-0.5*ewald_gamma))*exp(-ewald_gamma*ewc7))/ewc7

    der_grav_integral = 0.5D0*pi*(ERFC_SCALED(ewc2*(ewc7*ewald_dzeta+0.5*ewald_gamma)) &
    &             *exp(-ewald_gamma**2/(4.0D0*ewald_dzeta)-ewald_dzeta*ewc7*ewc7) - &
    &             ERFC(ewc2*(ewc7*ewald_dzeta-0.5*ewald_gamma))*exp(-ewald_gamma*ewc7))

    der2_grav_integral = (pi**2)*ewc2*Linv*(exp(-ewald_dzeta*ewc7*ewc7)*exp(-ewc6*ewald_gamma*ewald_gamma)* &
    &    (2.0D0/sqrt(pi) - (sqrt(ewald_dzeta)*ewc7)*ERFC_SCALED(ewc2*(ewc7*ewald_dzeta+0.5*ewald_gamma))) - &
    &    (sqrt(ewald_dzeta)*ewc7)*exp(-ewald_gamma*ewc7)*ERFC(ewc2*(ewc7*ewald_dzeta-0.5*ewald_gamma)))

! potential
    grv_bhAccLong2P1I(0) = grav_integral*ewald_betac*Linv2
! acceleration
    grv_bhAccLong2P1I(1) = 4*hi*ratio_pinv*(Linv**2) * ewald_betas*grav_integral
    grv_bhAccLong2P1I(2) = 4*hjr*ratio_pinv*(Linv**2) * ewald_betas*grav_integral
    grv_bhAccLong2P1I(3) = - 2*ratio_pinv*(Linv**2)*ewald_betac*der_grav_integral
! partial derivatives of acceleration
    grv_bhAccLong2P1I(4) =  8*pi*ratio_pinv*(Linv**3)*(hi**2)*ewald_betac*grav_integral
    grv_bhAccLong2P1I(5) = 8*pi*(ratio_pinv**2)*(Linv**3)*hi*hjr*ewald_betac*grav_integral
    grv_bhAccLong2P1I(6) = 4*pi*ratio_pinv*(Linv**3)*hi*ewald_betas*der_grav_integral
    grv_bhAccLong2P1I(8) =  8*pi*(ratio_pinv**3)*(Linv**3)*(hjr**2)*ewald_betac*grav_integral
    grv_bhAccLong2P1I(9) =  4*pi*(ratio_pinv**2)*(Linv**3)*hjr*ewald_betas*der_grav_integral
    grv_bhAccLong2P1I(12) =  2*ratio_pinv*(Linv**2)*ewald_betac*der2_grav_integral
  endif

end


