!!****if* source/Simulation/SimulationMain/unitTest/ProtonImaging/Validation1/Simulation_data
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

#include "ProtonImaging.h"
#include "constants.h"

  character (len = MAX_STRING_LENGTH), save :: sim_baseName

  logical, save :: sim_printBlockVariables
  logical, save :: sim_clockwiseB

  integer, save :: sim_globalComm
  integer, save :: sim_globalMe
  integer, save :: sim_globalNumProcs
  integer, save :: sim_nCircle
  integer, save :: sim_refinementLevel

  real,    save :: sim_cellSizeX
  real,    save :: sim_cellSizeY
  real,    save :: sim_cellSizeZ
  real,    save :: sim_magneticFluxDensity   ! in g/(esu * s)
  real,    save :: sim_protonCharge
  real,    save :: sim_protonMass
  real,    save :: sim_speedOfLight
  real,    save :: sim_xCenter
  real,    save :: sim_zCenter

  real,    save :: sim_xCircle (1:BEAM_GRIDARRAYSIZE)
  real,    save :: sim_yCircle (1:BEAM_GRIDARRAYSIZE)

  real, allocatable :: sim_screenX (:)
  real, allocatable :: sim_screenY (:)

end module Simulation_data
