!!****if* source/physics/materialProperties/MagneticResistivity/MagneticResistivityMain/SpitzerHighZ/MagneticResistivity_init
!!
!! NAME
!!
!!  MagneticResistivity_init
!!
!! SYNOPSIS
!!
!!  MagneticResistivity_init()
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!
!!
!!
!!***
subroutine MagneticResistivity_init
  use MagneticResistivity_data, ONLY: mag_useMagneticResistivity
  use MagneticResistivity_data, ONLY: res_mele
  use MagneticResistivity_data, ONLY: res_qele
  use MagneticResistivity_data, ONLY: res_navo
  use MagneticResistivity_data, ONLY: res_hbar
  use MagneticResistivity_data, ONLY: res_boltz
  use MagneticResistivity_data, ONLY: res_speedlt
  use MagneticResistivity_data, ONLY: res_ieTimeCoef
  use MagneticResistivity_data, ONLY: res_mUnit
  use MagneticResistivity_data, ONLY: res_coef
  use MagneticResistivity_data, ONLY: res_maxRes

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  use PhysicalConstants_interface, ONLY : PhysicalConstants_get

  implicit none

#include "Flash.h"
#include "constants.h"

  call RuntimeParameters_get("useMagneticResistivity", mag_useMagneticResistivity)
  call RuntimeParameters_get("res_ieTimeCoef", res_ieTimeCoef)
  call RuntimeParameters_get("UnitSystem", res_mUnit)
  call RuntimeParameters_get("res_coef", res_coef)
  call RuntimeParameters_get("res_maxRes", res_maxRes)
  if (res_maxRes < 0.0) then
     res_maxRes = HUGE(res_maxRes) ! Disable ceiling if RP res_maxRes is given as a negative number
  end if

  ! Set physical constants:
  call PhysicalConstants_get("electron mass",res_mele)
  call PhysicalConstants_get("electron charge",res_qele)
  call PhysicalConstants_get("speed of light",res_speedlt)
  call PhysicalConstants_get("Avogadro", res_navo)
  call PhysicalConstants_get("Boltzmann",res_boltz)
  call PhysicalConstants_get("Planck",res_hbar)
  res_hbar = res_hbar/(2.0*PI)

end subroutine MagneticResistivity_init

