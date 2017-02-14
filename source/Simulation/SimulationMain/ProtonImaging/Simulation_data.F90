!!****if* source/Simulation/SimulationMain/ProtonImaging/Simulation_data
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
!!  Stores the local data for the proton imaging unit test.
!!  
!!***

Module Simulation_data

  implicit none

#include "constants.h"

  character (len = MAX_STRING_LENGTH), save :: sim_baseName

  logical, save :: sim_printBlockVariables

  integer, save :: sim_globalComm
  integer, save :: sim_globalMe
  integer, save :: sim_globalNumProcs

  real,    save :: sim_cellSizeX
  real,    save :: sim_cellSizeY
  real,    save :: sim_cellSizeZ

end module Simulation_data
