!!****if* source/physics/materialProperties/MagneticResistivity/MagneticResistivityMain/SpitzerHighZ/res_ieEquilTime
!!
!! NAME
!!
!!  res_ieEquilTime
!!
!! SYNOPSIS
!!
!!  call res_ieEquilTime(real, intent(in )  :: zbar,
!!                      real, intent(in )  :: abar,
!!                      real, intent(in )  :: tele,
!!                      real, intent(in )  :: tion,
!!                      real, intent(in )  :: nion,
!!                      real(out) :: eqtime)
!!
!! DESCRIPTION
!!
!! Compute the Coulomb Logarithm using the form presented in:
!!     H. Brysk, P.M. Campbell, P. Hammerling, 
!!     "Thermal Conduction in Laser Fusion", 
!!     Plasma Physics, vol. 17, pp. 473, (1974). 
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

subroutine res_ieEquilTime(zbar, abar, tele, tion, nion, eqtime)
  use MagneticResistivity_data,ONLY: res_ieTimeCoef, res_boltz, res_qele, &
       res_navo, res_mele, res_hbar
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

  call loglambda(tele, zbar*nion, zbar, ll)
  mion = abar / res_navo

  eqtime = res_ieTimeCoef * &
       3.0 * res_boltz**1.5 / (8.0 * sqrt(2*PI) * res_qele**4) * &
       (mion * tele + res_mele * tion)**1.5 / &
       ( sqrt(mion*res_mele) * zbar**2 * nion * ll)

contains

  subroutine loglambda(tele, nele, zbar, ll)
    implicit none

    ! This subroutine computes the Coulomb logarithm. The formula used
    ! comes from Atzeni.
    
    real, intent(in)  :: tele ! electron temperature [K]
    real, intent(in)  :: nele ! electron number density [cm^-3]
    real, intent(in)  :: zbar ! the average ionization [unitless]
    real, intent(out) :: ll   ! the coulomb logarithm [unitless]

    real :: bmax, bmin, bmin_classic, bmin_quantum
    real, parameter :: ll_floor = 1.0
    
    bmax = sqrt(res_boltz * tele / (4*PI * res_qele**2 * nele))
    
    bmin_classic = zbar * res_qele**2 / (3*res_boltz*tele)
    bmin_quantum = res_hbar / (2*sqrt(3*res_boltz*tele*res_mele))
    bmin = max(bmin_classic, bmin_quantum)

    ll = max(log(bmax/bmin), ll_floor)

  end subroutine loglambda


end subroutine res_ieEquilTime
