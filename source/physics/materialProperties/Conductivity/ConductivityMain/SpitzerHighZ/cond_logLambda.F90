!!****if* source/physics/materialProperties/Conductivity/ConductivityMain/SpitzerHighZ/cond_logLambda
!!
!! NAME
!!
!!  cond_logLambda
!!
!! SYNOPSIS
!!
!!  call cond_logLambda(real(IN)  :: tele,
!!                      real(IN)  :: nele,
!!                      real(IN)  :: zbar,
!!                      real(OUT) :: ll)
!!
!! DESCRIPTION
!!
!! Computes the log of Lambda (Coulomb logarithm) for Spitzer electron conductivity.
!! USed by SpitzerHighZ Conductivity, see there for further comments.
!!
!! ARGUMENTS
!!
!!   tele      :   electron temperature
!!   nele      :   electron number density
!!   zbar      :   Zbar
!!   ll        :   log of Lambda
!!
!! SEE ALSO
!!
!!  Conductivity 
!!***

#include "constants.h"

  subroutine cond_logLambda(tele, nele, zbar, ll)
    use Conductivity_data, ONLY: cond_mele, cond_boltz, cond_hbar, cond_qele
    implicit none

    ! This subroutine computes the Coulomb logarithm. The formula used
    ! comes from Atzeni.

    real, intent(in)  :: tele ! electron temperature [K]
    real, intent(in)  :: nele ! electron number density [cm^-3]
    real, intent(in)  :: zbar ! the average ionization [unitless]
    real, intent(out) :: ll   ! the coulomb logarithm [unitless]

    real :: bmax, bmin, bmin_classic, bmin_quantum
    real, parameter :: ll_floor = 1.0

    bmax = sqrt(cond_boltz * tele / (4*PI * cond_qele**2 * nele))

    bmin_classic = zbar * cond_qele**2 / (3*cond_boltz*tele)
    bmin_quantum = cond_hbar / (2*sqrt(3*cond_boltz*tele*cond_mele))
    bmin = max(bmin_classic, bmin_quantum)

    ll = max(log(bmax/bmin), ll_floor)

  end subroutine cond_logLambda

