!!****if* source/physics/materialProperties/MagneticResistivity/MagneticResistivityMain/SpitzerHighZ/MagneticResistivity_fullState
!!
!! NAME
!!  MagneticResistivity_fullState
!!
!! SYNOPSIS
!!  call MagneticResistivity_fullState(real(in)    :: solnVec(NUNK_VARS),
!!                                     real(out)   :: resPar,
!!                            OPTIONAL,real(out)   :: resPerp)
!!
!! DESCRIPTION
!!
!! Computes the Spitzer electron Magnetic Resistivity for all materials,
!! including those with Z > 1. The specific equations used here all
!! come from "The Physics of Inertial Fusion" by Atzeni.
!!
!!  Returns Magnetic Resistivity, parallel and perpendicular components.
!!
!! ARGUMENTS
!!
!!   solnVec  :   solution state, a vector from UNK with all variables
!!   resPar   :   parallel component of Magnetic Resistivity
!!   resPerp :    perpendicular component of Magnetic Resistivity
!!
!!***

#include "Flash.h"
#include "constants.h"  


subroutine MagneticResistivity_fullState(solnVec,resPar, resPerp)
  use MagneticResistivity_interface, ONLY: MagneticResistivity
  use MagneticResistivity_data, ONLY: mag_useMagneticResistivity, &
       res_mele, res_qele, res_navo, res_speedlt
  use MagneticResistivity_data, ONLY: res_mUnit
  use MagneticResistivity_data, ONLY: res_coef
  use MagneticResistivity_data, ONLY: res_maxRes
  use Eos_interface, ONLY: Eos_getAbarZbar

  implicit none

  real, intent(in)  :: solnVec(NUNK_VARS)
  real, intent(out) :: resPar
  real, OPTIONAL, intent(out) :: resPerp

  real :: dens
  real :: tele
  real :: tion
  real :: nele
  real :: nion
  real :: abar
  real :: zbar
  real :: eqtime
  real :: resPerpLoc

  call Eos_getAbarZbar(solnVec=solnVec,abar=abar,zbar=zbar)

  dens = solnVec(DENS_VAR)
  nion = dens * res_navo / abar
  nele = zbar * nion

#ifdef FLASH_3T
  tele = solnVec(TELE_VAR)
  tion = solnVec(TION_VAR)
#else
  tele = solnVec(TEMP_VAR)
  tion = solnVec(TEMP_VAR)
#endif

  call res_ieEquilTime(zbar, abar, tele, tion, nion, eqtime)

  resPerpLoc = res_coef * res_mele / (res_qele**2 * nele * eqtime)
  if (present(resPerp)) resPerp = resPerpLoc
  ! This formula is only valid when the magnetic field is strong (and
  ! only for hydrogen):
  resPar = resPerpLoc/1.96

  !! Scale resistivity depending on units
  if (res_mUnit == "SI" .or. res_mUnit == "si" ) then
     resPerpLoc = resPerpLoc*1.e7/(4.*PI)
     resPar  = resPar *1.e7/(4.*PI)
  elseif (res_mUnit == "CGS" .or. res_mUnit == "cgs" ) then
     resPerpLoc = resPerpLoc*res_speedlt**2/(4.*PI)
     resPar  = resPar *res_speedlt**2/(4.*PI)
  else !no unit
     ! Do nothing
     resPerpLoc = resPerpLoc*res_speedlt**2/(4.*PI)
     resPar  = resPar *res_speedlt**2/(4.*PI)
  end if

  if (present(resPerp)) resPerp = resPerpLoc
  !! cap resPar and resPerp to avoid unphysically cold regions
  if (present(resPerp)) resPerp = min(resPerp, res_maxRes)
  resPar = min(resPar, res_maxRes)
 

end subroutine MagneticResistivity_fullState
