!!****if* source/physics/sourceTerms/Heatexchange/HeatexchangeMain/LeeMore/Heatexchange_data
!!
!! NAME
!!
!!  Heatexchange_data
!!
!! SYNOPSIS
!!  use Heatexchange_data
!!
!! DESCRIPTION
!!
!!  Defines and stores the local data for Unit Heatexchange. All the variables
!!  defined here are initialized in Heatexchange_init() by calling 
!!  RuntimeParameters_get subroutine. These data variables are for 
!!  Unit Scope ->  Heatexchange, and can be used by all Heatexchnage subroutines.
!!  
!!***


Module Heatexchange_data

  implicit none

#include "constants.h"

  logical, save :: hx_useHeatexchange
  logical, save :: hx_restart
  integer, save :: hx_meshMe
  
  ! Coefficient for ion/electron coupling time
  real,    save :: hx_ieTimeCoef

  ! Physical constants:
  real, save :: hx_mele  ! Electron mass (grams)
  real, save :: hx_boltz ! Boltzmann constant (ergs/K)
  real, save :: hx_hbar  ! Planck's constant over 2 PI (erg*s)
  real, save :: hx_qele  ! Proton charge (esu)
  real, save :: hx_navo  ! Avogadros number

  ! These variables are used to compute the hx DT:
  real, save :: hx_smallE
  real, save :: hx_dtFactor
  real, save :: hx_relTol

end Module Heatexchange_data
