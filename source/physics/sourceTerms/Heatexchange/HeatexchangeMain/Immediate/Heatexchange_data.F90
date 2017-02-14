!!****if* source/physics/sourceTerms/Heatexchange/HeatexchangeMain/Immediate/Heatexchange_data
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
!!  Defines and stores the local data for Unit Heatexchange. The variables
!!  defined here are initialized in Heatexchange_init(), mostly by calling the
!!  RuntimeParameters_get subroutine. These data variables are for 
!!  Unit Scope -> Heatexchange, and are specific to the "Immediate" implementation.
!!  
!!***


Module Heatexchange_data

  implicit none

  integer, save :: hx_globalMe
  logical, save :: hx_useHeatexchange, hx_restart
  logical, save :: hx_applyToRadiation

end Module Heatexchange_data
