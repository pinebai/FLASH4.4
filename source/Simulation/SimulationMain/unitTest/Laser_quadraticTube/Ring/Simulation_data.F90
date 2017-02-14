!!****if* source/Simulation/SimulationMain/unitTest/Laser_quadraticTube/Ring/Simulation_data
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
  real,    save :: sim_maxDeltaP
  real,    save :: sim_maxDeltaRsqr
  real,    save :: sim_nc
  real,    save :: sim_nuw
  real,    save :: sim_nw
  real,    save :: sim_powerDecayFactor
  real,    save :: sim_R
  real,    save :: sim_Tw
  real,    save :: sim_xc
  real,    save :: sim_xw
  real,    save :: sim_yfocal
  real,    save :: sim_Z
  real,    save :: sim_zc
  real,    save :: sim_zw

  integer, parameter :: sim_totalNumberOfBoxes  =  10000
  integer, parameter :: sim_refinementLevelMax  =  6
  real,    parameter :: sim_keV2Kelvin          =  11604505.

  integer, save :: sim_powerBox  (1:sim_totalNumberOfBoxes)
  integer, save :: sim_radialBox (1:sim_totalNumberOfBoxes)
  real,    save :: sim_boxPower  (1:sim_totalNumberOfBoxes)
  real,    save :: sim_boxRadius (1:sim_totalNumberOfBoxes)

end module Simulation_data
