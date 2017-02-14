!!****if* source/Simulation/SimulationMain/magnetoHD/Torus/Simulation_data
!!
!! NAME
!!
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  Simulation_data()
!!
!!  Stores the local data for Simulation setup: 2D Torus MHD problem
!!
!!  
!!
!!***

module Simulation_data

  implicit none
#include "constants.h"

  !! *** Runtime Parameters *** !!
  real, save    :: sim_gamma
  real, save    :: sim_smallP
  real, save    :: sim_xMin, sim_xMax, sim_yMin, sim_yMax, sim_zMin, sim_zMax
  real, save    :: sim_rMin, sim_rMax, sim_denMax, sim_denCut, sim_beta, sim_dCon
  real, save    :: sim_tCon, sim_r0, sim_rSphere

  logical, save :: sim_gCell, sim_killdivb

  integer, save :: sim_meshMe
end module Simulation_data

