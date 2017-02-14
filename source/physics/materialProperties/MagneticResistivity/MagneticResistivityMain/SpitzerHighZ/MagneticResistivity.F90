!!****if* source/physics/materialProperties/MagneticResistivity/MagneticResistivityMain/SpitzerHighZ/MagneticResistivity
!!
!! NAME
!!  MagneticResistivity
!!
!! SYNOPSIS
!!  call MagneticResistivity  (real(in)  :: dens,
!!                            real(in)   :: temp,
!!                            real (in)  :: xn(NSPECIES) ,
!!                            real(out)  :: magResist)
!!
!! DESCRIPTION
!!
!! Computes the Spitzer electron Magnetic Resistivity for all materials,
!! including those with Z > 1. The specific equations used here all
!! come from "The Physics of Inertial Fusion" by Atzeni.
!!
!!  Returns Magnetic Resistivity resPar in magResist.
!!
!! ARGUMENTS
!!
!! ARGUMENTS
!!  temp      - Plasma temperature
!!  dens      - Plasma density
!!  xn        - Species
!!  magResist - Magnetic resistivity
!!
!!***

#include "Flash.h"
#include "constants.h"  

subroutine MagneticResistivity(temp,dens,xn,magResist)
  use MagneticResistivity_data, ONLY: mag_useMagneticResistivity, &
       res_mele, res_qele, res_navo, res_speedlt
  use MagneticResistivity_data, ONLY: res_mUnit
  use MagneticResistivity_data, ONLY: res_coef
  use MagneticResistivity_data, ONLY: res_maxRes
  use Eos_interface, ONLY: Eos_getAbarZbar

  implicit none

  real, intent(IN) :: temp, dens
  real, intent(IN), dimension(NSPECIES) :: xn
  real, intent(OUT):: magResist

  real :: resPar
  real :: tele
  real :: tion
  real :: nele
  real :: nion
  real :: abar
  real :: zbar
  real :: eqtime
  real :: resPerpLoc

  call Eos_getAbarZbar(massFrac=xn,abar=abar,zbar=zbar)

  nion = dens * res_navo / abar
  nele = zbar * nion

  tele = temp
  tion = temp

  call res_ieEquilTime(zbar, abar, tele, tion, nion, eqtime)

  resPerpLoc = res_coef * res_mele / (res_qele**2 * nele * eqtime)

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
  !! cap resPar and resPerp to avoid unphysically cold regions
  resPar = min(resPar, res_maxRes)
  
  magResist = resPar
  
end subroutine MagneticResistivity
