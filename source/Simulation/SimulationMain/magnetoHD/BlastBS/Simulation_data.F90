!!****if* source/Simulation/SimulationMain/magnetoHD/BlastBS/Simulation_data
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
!!  Stores the local data for Simulation setup: BlastBS
!!  
!!
!!
!!***

module Simulation_data

  implicit none
#include "constants.h"
#include "Flash.h"

  integer, save :: sim_meshMe

  !! *** Runtime Parameters *** !!
  real, save :: sim_gamma
  real, save :: sim_smallX, sim_smallRho, sim_smallP
  real, save :: sim_xMin, sim_xMax, sim_yMin, sim_yMax, sim_zMin, sim_zMax
  real, save :: sim_xCtr, sim_yCtr, sim_zCtr
  real, save :: sim_Radius
  logical, save :: sim_gCell

#if defined(MAGX_VAR) && defined(MAGY_VAR) && defined(MAGZ_VAR)
  real, save :: sim_Bx0
  logical, save :: sim_killdivb
#endif

end module Simulation_data
