!!****if* source/Simulation/SimulationMain/unitTest/Laser_quadraticTube/testII/Simulation_data
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
!!  Stores the local data for the laser quadratic tube unit test.
!!  
!!***

Module Simulation_data

  implicit none

#include "constants.h"

  character (len = MAX_STRING_LENGTH), save :: sim_baseName

  logical, save :: sim_printBlockVariables

  integer, save :: sim_geometry
  integer, save :: sim_globalComm
  integer, save :: sim_globalMe
  integer, save :: sim_globalNumProcs
  integer, save :: sim_numRaysAnalyzed
  integer, save :: sim_numRaysLaunched
  integer, save :: sim_refinementLevel

  real,    save :: sim_A
  real,    save :: sim_beamFrequency
  real,    save :: sim_beamTargetRadius
  real,    save :: sim_beamWavelength
  real,    save :: sim_focalPointRadius
  real,    save :: sim_maxDeltaPZ
  real,    save :: sim_maxDeltaRsqr
  real,    save :: sim_nc
  real,    save :: sim_nuw
  real,    save :: sim_nw
  real,    save :: sim_percentFocus
  real,    save :: sim_powerDecayFactor
  real,    save :: sim_powerPartition
  real,    save :: sim_tcross
  real,    save :: sim_Tw
  real,    save :: sim_xc
  real,    save :: sim_xw
  real,    save :: sim_yfocal
  real,    save :: sim_Z
  real,    save :: sim_zc
  real,    save :: sim_zw

  integer, parameter :: GRID_3DCARTESIAN    =  1, &
                        GRID_2DCYLINDRICAL  =  2, &
                        GRID_2DCARTESIAN    =  3

  integer, parameter :: sim_refinementLevelMax  =  6
  real,    parameter :: sim_keV2Kelvin          =  11604505.

  real,    save :: sim_maxDeltaRsqrErrorBar (1:sim_refinementLevelMax)
  real,    save :: sim_maxDeltaPZErrorBar   (1:sim_refinementLevelMax)

end module Simulation_data
