!!****if* source/physics/materialProperties/Conductivity/ConductivityMain/SpitzerHighZ/Conductivity_data
!!
!! NAME
!!
!!  Conductivity_data
!!
!! SYNOPSIS
!!
!!  use Conductivity_data
!!
!! DESCRIPTION
!!
!!  Stores the local data for the Spitzer Conductivity implementation.
!!
!!
!!***

Module Conductivity_data

  implicit none

  logical, save::  cond_useConductivity
  
  ! Storage for some physical constants:
  real, save :: cond_mele  ! Electron mass (grams)
  real, save :: cond_boltz ! Boltzmann constant (ergs/K)
  real, save :: cond_hbar  ! Planck's constant over 2 PI (erg*s)
  real, save :: cond_qele  ! Proton charge (esu)
  real, save :: cond_navo  ! Avogadros number
  
end Module Conductivity_data
