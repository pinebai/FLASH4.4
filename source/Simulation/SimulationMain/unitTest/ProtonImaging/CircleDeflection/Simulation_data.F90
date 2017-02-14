!!****if* source/Simulation/SimulationMain/unitTest/ProtonImaging/CircleDeflection/Simulation_data
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
  logical, save :: sim_clockwiseB

  integer, save :: sim_globalComm
  integer, save :: sim_globalMe
  integer, save :: sim_globalNumProcs
  integer, save :: sim_refinementLevel

!  integer, parameter :: sim_nBeamPoints = 6825  ! square -> 4 x 1000 , cross -> 2 x 1412 , center -> 1
  integer, parameter :: sim_nBeamPoints = 10000  ! for a circle

  real,    save :: sim_cellSizeX
  real,    save :: sim_cellSizeY
  real,    save :: sim_cellSizeZ
  real,    save :: sim_magneticFluxDensity   ! in g/(esu * s)
  real,    save :: sim_protonCharge
  real,    save :: sim_protonMass
  real,    save :: sim_speedOfLight
  real,    save :: sim_xCenter
  real,    save :: sim_zCenter

  real,    save :: sim_beamPointsX (1:sim_nBeamPoints)
  real,    save :: sim_beamPointsY (1:sim_nBeamPoints)

  real, allocatable :: sim_screenX (:)
  real, allocatable :: sim_screenY (:)

end module Simulation_data
