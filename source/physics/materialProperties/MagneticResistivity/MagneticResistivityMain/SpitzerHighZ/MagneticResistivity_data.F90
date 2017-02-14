!!****if* source/physics/materialProperties/MagneticResistivity/MagneticResistivityMain/SpitzerHighZ/MagneticResistivity_data
!!
!! NAME
!!
!!  MagneticResistivity_data
!!
!! SYNOPSIS
!!
!!  use MagneticResistivity_data
!!
!! DESCRIPTION
!!
!!  Stores the local data for the Spitzer MagneticResistivity implementation.
!!
!!
!!***

Module MagneticResistivity_data

  implicit none

  logical, save::  mag_useMagneticResistivity
  
  ! Storage for some physical constants:
  real, save :: res_mele       ! Electron mass (grams)
  real, save :: res_qele       ! Proton charge (esu)
  real, save :: res_navo       ! Avogadros number
  real, save :: res_speedlt    ! Speed of light (cm/s)
  real, save :: res_boltz      ! Boltzmann's constant (ergs/K)
  real, save :: res_hbar       ! Planck's constant / 2 / PI (ergs * s)
  real, save :: res_ieTimeCoef ! Ion/electron equilibration time coefficient
  real, save :: res_coef       ! Resistivity scaling coefficient
  real, save :: res_maxRes     ! Ceiling value

  character(4), save :: res_mUnit ! Unit system

end Module MagneticResistivity_data
