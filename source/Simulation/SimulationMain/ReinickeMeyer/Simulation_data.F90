!!****if* source/Simulation/SimulationMain/ReinickeMeyer/Simulation_data
!!
!! NAME
!!
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  use Simulation_data 
!!
!!  DESCRIPTION
!!
!!  Stores the local data for Simulation setup
!!  
!!***

module Simulation_data

  implicit none
#include "constants.h"

  !! *** Runtime Parameters *** !!
  real, save :: sim_singleSpeciesA
  real, save :: sim_DensityExponent
  real, save :: sim_TemperatureExponent 
  real, save :: sim_K0
  real, save :: sim_gamma
  real, save :: sim_tInitial
  real, save :: sim_smlrho
  real, save :: sim_smallt

  ! Initial thermal front position [cm]:
  real, save :: sim_rfInit
  
  !! *** Physical Constants *** !!

  ! Avogadros number:
  real, save :: sim_avo

  ! Boltzmann's constant:
  real, save :: sim_boltz

  ! Geometry:
  integer :: sim_geometry

end module Simulation_data
