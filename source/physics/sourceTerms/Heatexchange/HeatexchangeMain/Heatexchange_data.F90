!!****if* source/physics/sourceTerms/Heatexchange/HeatexchangeMain/Heatexchange_data
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
!!  Defines and stores the local data for Unit Heatexchange. Many of the variables
!!  defined here are initialized in Heatexchange_init() by calling the
!!  RuntimeParameters_get subroutine. These data variables are for 
!!  Unit Scope ->  Heatexchange, and can be used by various Heatexchange subroutines.
!!  
!!***


Module Heatexchange_data

  implicit none

#include "constants.h"

  integer, save :: hx_meshMe
  logical, save :: hx_useHeatexchange, hx_restart

  integer, save :: hx_logLevel = 700
  
  integer, save :: hx_geometry
  real,    save :: hx_smallE
  real,    save :: hx_c12, hx_c13, hx_c23
  real,    save :: hx_coulombLog
  real,    save :: hx_singleSpeciesA, hx_singleSpeciesZ
  real,    save :: hx_Avogadro, hx_kBoltzmann, hx_eleCharge, hx_eMassInUAmu

  real, save :: hx_dtFactor
  real, save :: hx_relTol

  ! Flag to indicate whether the RadTrans unit is in use 
  logical, save :: hx_useRadTrans 

end Module Heatexchange_data
