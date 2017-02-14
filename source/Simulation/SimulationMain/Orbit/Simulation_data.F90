!!****if* source/Simulation/SimulationMain/Orbit/Simulation_data
!!
!! NAME
!!
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  use Simulation_data
!!
!! DESCRIPTION
!!
!!  Stores the local data for Simulation setup:  Orbit
!!  
!!
!!***

module Simulation_data
  
  implicit none
  
#include "constants.h"
  
  !! *** Runtime Parameters *** !!
  
  real, save    :: sim_xMin, sim_xMax, sim_yMin, sim_yMax, sim_zMin, sim_zMax
  real, save    :: sim_Sep, sim_ptMass
  logical, save :: sim_extField
  
  !! *** Variables pertaining to this Simulation *** !!

  integer, save      :: sim_nPtot
  integer, save      :: sim_NumPEs, ptagmult
  ! Small differences in x, y, z directions
  real, save         :: sim_EPSX, sim_EPSY, sim_EPSZ
  real, save         :: sim_xPosns(2)
  real, save         :: sim_yPosns(2)
  real, save         :: sim_zPosns(2)
  real, parameter    :: EPS = 1.0E-6
  real, save         :: sim_Newton

integer, save :: sim_meshMe
end module Simulation_data
