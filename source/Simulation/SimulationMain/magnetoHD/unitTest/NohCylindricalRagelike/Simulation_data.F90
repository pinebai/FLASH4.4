!!****if* source/Simulation/SimulationMain/magnetoHD/unitTest/NohCylindricalRagelike/Simulation_data
!!
!! NAME
!!
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  Simulation_data()
!!
!!  Stores the local data for Simulation setup: Cylindrical Noh
!!  
!!  Reference: Velikovich et al. Phys. Plasmas 19 (2012), 012707
!!
!!***

module Simulation_data

  implicit none
#include "constants.h"

  !! *** Runtime Parameters *** !!
  real, save    :: sim_gamma
  real, save    :: sim_smallP
  real, save    :: sim_xMin, sim_xMax, sim_yMin, sim_yMax, sim_zMin, sim_zMax
  real, save    :: sim_UnitDensity
  real, save    :: sim_UnitVelocity
  real, save    :: sim_UnitLength

  logical, save :: sim_gCell, sim_killdivb

integer, save :: sim_meshMe
end module Simulation_data
