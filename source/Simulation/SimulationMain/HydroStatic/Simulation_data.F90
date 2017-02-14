!!****if* source/Simulation/SimulationMain/HydroStatic/Simulation_data
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
!!  Stores the local data for Simulation setup: Hydrostatic
!!  
!!
!!
!!***
module Simulation_data
  implicit none
#include "constants.h"
  real, save :: sim_gamma, sim_smallX
  integer, save :: sim_gravVector(MDIM)
  real, save :: sim_GravConst
  integer, save :: sim_gravDirec

  real, save :: sim_xyzRef, sim_presRef, sim_densRef, sim_tempRef
  real, save :: sim_molarMass, sim_gasconstant
integer, save :: sim_meshMe
end module Simulation_data
