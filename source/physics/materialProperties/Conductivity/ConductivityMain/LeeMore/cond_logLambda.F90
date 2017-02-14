!!****if* source/physics/materialProperties/Conductivity/ConductivityMain/LeeMore/cond_logLambda
!!
!! NAME
!!
!!  cond_logLambda
!!
!! SYNOPSIS
!!
!!  call cond_logLambda(real(IN)  :: tele,
!!                      real(IN)  :: nele,
!!                      real(IN)  :: tion,
!!                      real(IN)  :: zbar,
!!                      real(OUT) :: ll)
!!
!! DESCRIPTION
!!
!! Computes the log of Lambda (Coulomb logarithm) for the Lee & More
!! (Phys. Fluids 27, 1273 1984) electron conductivity model. 
!!
!! ARGUMENTS
!!
!!   tele      :   electron temperature
!!   nele      :   electron number density
!!   tion      :   ion temperature
!!   zbar      :   Zbar
!!   ll        :   log of Lambda
!!
!! SEE ALSO
!!
!!  Conductivity 
!!***

#include "constants.h"

  subroutine cond_logLambda(tele, nele, tion, zbar, ll)
    use Conductivity_data, ONLY: cond_mele, cond_boltz, cond_hbar, cond_qele
    implicit none

    ! This subroutine computes the Coulomb logarithm. The formulae used
    ! come from Lee & More 1984.

    real, intent(in)  :: tele ! electron temperature [K]
    real, intent(in)  :: nele ! electron number density [cm^-3]
    real, intent(in)  :: tion ! ion temperature [K]
    real, intent(in)  :: zbar ! the average ionization [unitless]
    real, intent(out) :: ll   ! the coulomb logarithm [unitless]

    real :: bmax, bmin, bmin_classic, bmin_quantum
    real :: debye_len, nion, R0, Tf
    real, parameter :: ll_floor = 2.0

    ! Calculate interatomic distance
    nion = nele/zbar
    R0 = (4*PI*nion/3)**(-1./3)

    !Fermi temperature in Kelvin
    Tf = cond_hbar**2 / (2*cond_mele) * (3 * PI**2 * nele)**(2./3) / cond_boltz; 

    ! LM84 Eq. 19
    debye_len = 1/sqrt((4*PI * cond_qele**2 * nele)/(cond_boltz*sqrt(tele**2 + Tf**2)) + &
        4*PI*nion *zbar**2 *cond_qele**2 /(cond_boltz*tion) )

    bmax = max(debye_len,R0)

    bmin_classic = zbar * cond_qele**2 / (3*cond_boltz*tele)

    ! This correctly uses hbar on the top instead of just h (as in LM84)
    ! as pointed out by Skupsky, S. Phys. Rev. A. 36, 5701 (1987)
    bmin_quantum = cond_hbar / (2*sqrt(3*cond_boltz*tele*cond_mele))

    bmin = max(bmin_classic, bmin_quantum)
    
    ll = max(0.5*log(1+(bmax/bmin)**2), ll_floor)


  end subroutine cond_logLambda

