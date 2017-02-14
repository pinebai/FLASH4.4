!!****if* source/physics/Hydro/HydroMain/split/RHD/hy_rhd_pressureFunc
!!
!! NAME
!!
!!   hy_rhd_pressureFunc
!!
!! SYNOPSIS
!!   hy_rhd_pressureFunc( real    (IN) :: u(:),
!!                        real    (IN) :: x(:),
!!                        real   (OUT) :: fp(:),
!!                        real   (OUT) :: dfdp(:),
!!                        integer(OUT) :: iflag)
!! 
!!
!!
!!    Define the [ f(p) - p]  function, and its derivative,
!!    whose zero gives the pressure.
!!    X is the independent variable (pressure)
!!
!!   u : vector of conservative quantities
!!
!! ARGUMENTS
!! 
!!  u     -      a vector of primitive quantities
!!  x     -      independent variable (pressure)
!!  fp    -      the storage for the pressure func
!!  dfdp  -      storage for the derivative of fp
!! iflag  -      used for debugging 
!!
!!***
subroutine hy_rhd_pressureFunc(u, x, fp, dfdp, iflag)

  use Hydro_data, ONLY : hy_gamma

  implicit none

#include "Flash.h"
#include "RHD.h"  
  
  integer, INTENT(out)                  :: iflag
  real, DIMENSION(NUNK_VARS),INTENT(in) :: u
  real, INTENT(in)                      :: x
  real, INTENT(out)                     :: fp, dfdp

  real  :: alpha, alpha2, S2
  real  :: lorentz, lorentz2
  real  :: tau, theta, h, scrh
  real  :: dh_dp, dh_dtau, gamma_r

  gamma_r  = hy_gamma/(hy_gamma - 1.d0)
  S2       = u(VELX_VAR)*u(VELX_VAR) &
           + u(VELY_VAR)*u(VELY_VAR) &
           + u(VELZ_VAR)*u(VELZ_VAR)
  alpha    = u(ENER_VAR) + x
  alpha2   = alpha*alpha
  lorentz2 = 1.e0 - S2/alpha2

  lorentz2 = max(lorentz2,1.e-9)

!  if (lorentz2 < 0.d0) then 
!    lorentz2 = 1.d-3
!    print *,' g < 0 in pressure fun'
!    iflag = 1
!    return
!  end  if

  iflag    = 0
  lorentz2 = 1.d0/lorentz2
  lorentz  = sqrt(lorentz2)

  tau   = lorentz/u(DENS_VAR)
  theta = x*tau

! get Enthalpy h and dh/dtau

  h       = 1.d0 + gamma_r*theta
  dh_dp   = gamma_r*tau
  dh_dtau = gamma_r*x

! get  f(P) and df(p)/dp

  fp   = u(DENS_VAR)*h*lorentz - u(ENER_VAR) - x;
  dfdp = u(DENS_VAR)*lorentz*dh_dp - S2*lorentz2*lorentz/(alpha2*alpha)*  &
          (lorentz*dh_dtau + u(DENS_VAR)*h) - 1.d0;
               
end subroutine hy_rhd_pressureFunc
