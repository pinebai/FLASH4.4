!!****if* source/physics/materialProperties/Conductivity/ConductivityMain/SpitzerHighZ/cond_logLambdaII
!!
!! NAME
!!
!!  cond_logLambdaII
!!
!! SYNOPSIS
!!
!!  call cond_logLambdaII(real(in) :: tion,
!!                        real(in) :: tele,
!!                        real(in) :: nele,
!!                        real(in) :: zbar,
!!                        real(out) :: ll)
!!
!! DESCRIPTION
!!
!! Computes log lambda
!! 
!! ARGUMENTS
!!
!!   tion : Ion temperature
!!
!!   tele : Electron temperature
!!
!!   nele : Electron density
!!
!!   zbar : zbar
!!
!!   ll : log lambda
!!
!!
!!
!!***

subroutine cond_logLambdaII(tion, tele, nele, zbar, ll)
#include "constants.h"

  use Conductivity_data, ONLY: cond_mele
  use Conductivity_data, ONLY: cond_boltz
  use Conductivity_data, ONLY: cond_hbar
  use Conductivity_data, ONLY: cond_qele
  implicit none

  ! This subroutine computes the Coulomb logarithm. The formula used
  ! comes from Atzeni.

  real, intent(in)  :: tion ! ion temperature [K]
  real, intent(in)  :: tele ! electron temperature [K]
  real, intent(in)  :: nele ! electron number density [cm^-3]
  real, intent(in)  :: zbar ! the average ionization [unitless]
  real, intent(out) :: ll   ! the coulomb logarithm [unitless]

  real :: bmax
  real :: bmin
  real, parameter :: ll_floor = 1.0

  bmax = sqrt(cond_boltz * tele / (4*PI * cond_qele**2 * nele))
  bmin = zbar**2 * cond_qele**2 / (3*cond_boltz*tion)

  ll = max(log(bmax/bmin), ll_floor)

end subroutine cond_logLambdaII
