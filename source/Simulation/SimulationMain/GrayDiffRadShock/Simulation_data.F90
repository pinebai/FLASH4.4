!!****if* source/Simulation/SimulationMain/GrayDiffRadShock/Simulation_data
!!
!! NAME
!!  Simulation_data
!!
!! SYNOPSIS
!!  Use Simulation_data
!!
!! DESCRIPTION
!!
!!  Store the simulation data
!!
!! 
!!***

module Simulation_data

  implicit none
  
#include "constants.h"

  !! *** Runtime Parameters *** !!
  real, save :: sim_rho
  real, save :: sim_temp
  real, save :: sim_M0
  real, save :: sim_P0
 
end module Simulation_data


