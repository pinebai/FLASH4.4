!!****if* source/physics/sourceTerms/Heatexchange/HeatexchangeMain/LeeMore/hx_ieEquilTime
!!
!! NAME
!!
!!  hx_ieEquilTime
!!
!! SYNOPSIS
!!
!!  call hx_ieEquilTime(real, intent(in )  :: zbar,
!!                      real, intent(in )  :: abar,
!!                      real, intent(in )  :: tele,
!!                      real, intent(in )  :: tion,
!!                      real, intent(in )  :: nion,
!!                      real(out) :: eqtime)
!!
!! DESCRIPTION
!!
!! Compute the Coulomb Logarithm using the form presented in:
!! Lee & More (Phys. Fluids 27, 1273, 1984)
!!
!! ARGUMENTS
!!
!!   zbar : Average ionization
!!   abar : Average atomic mass
!!   tele : electron temperature
!!   tion : ion temperature
!!   nion : electron number density
!!   eqtime : ion/electron equilibration time
!!
!!
!!***

subroutine hx_ieEquilTime(zbar, abar, tele, tion, nion, eqtime)
  use Heatexchange_data,ONLY: hx_ieTimeCoef, hx_boltz, hx_qele, &
       hx_navo, hx_mele, hx_hbar
  implicit none

#include "constants.h"

  ! Arguments:
  real, intent(in ) :: zbar
  real, intent(in ) :: abar
  real, intent(in ) :: tele
  real, intent(in ) :: tion
  real, intent(in ) :: nion
  real, intent(out) :: eqtime

  ! Local variables:
  real :: ll
  real :: mion

  call loglambda(tele, zbar*nion, tion, zbar, ll)
  mion = abar / hx_navo

  eqtime = hx_ieTimeCoef * &
       3.0 * hx_boltz**1.5 / (8.0 * sqrt(2*PI) * hx_qele**4) * &
       (mion * tele + hx_mele * tion)**1.5 / &
       ( sqrt(mion*hx_mele) * zbar**2 * nion * ll)

contains

  subroutine loglambda(tele, nele, tion, zbar, ll)
    implicit none

    ! This subroutine computes the Coulomb logarithm. 
    
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
    Tf = hx_hbar**2 / (2*hx_mele) * (3 * PI**2 * nele)**(2./3) / hx_boltz; 

    ! LM84 Eq. 19
    debye_len = 1/sqrt((4*PI * hx_qele**2 * nele)/(hx_boltz*sqrt(tele**2 + Tf**2)) + &
        4*PI*nion *zbar**2 *hx_qele**2 /(hx_boltz*tion) )

    bmax = max(debye_len,R0)
    
    bmin_classic = zbar * hx_qele**2 / (3*hx_boltz*tele)

    ! This correctly uses hbar on the top instead of just h (as in LM84)
    ! as pointed out by Skupsky, S. Phys. Rev. A. 36, 5701 (1987)
    bmin_quantum = hx_hbar / (2*sqrt(3*hx_boltz*tele*hx_mele))

    bmin = max(bmin_classic, bmin_quantum)

    ll = max(0.5*log(1+(bmax/bmin)**2), ll_floor)


  end subroutine loglambda


end subroutine hx_ieEquilTime
