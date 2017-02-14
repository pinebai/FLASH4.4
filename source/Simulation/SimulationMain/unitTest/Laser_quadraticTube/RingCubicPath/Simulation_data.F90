!!****if* source/Simulation/SimulationMain/unitTest/Laser_quadraticTube/RingCubicPath/Simulation_data
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
  integer, save :: sim_totalNumberOfBoxes

  real,    save :: sim_beamFrequency
  real,    save :: sim_beamTargetRadius
  real,    save :: sim_beamWavelength
  real,    save :: sim_focalPointRadius
  real,    save :: sim_maxDeltaRsqr
  real,    save :: sim_avgDeltaR
  real,    save :: sim_rmsDeltaR
  real,    save :: sim_sigmaDeltaR
  real,    save :: sim_avgDeltaX
  real,    save :: sim_rmsDeltaX
  real,    save :: sim_sigmaDeltaX
  real,    save :: sim_nc
  real,    save :: sim_nuw
  real,    save :: sim_nw
  real,    save :: sim_R
  real,    save :: sim_Tw
  real,    save :: sim_xc
  real,    save :: sim_xw
  real,    save :: sim_yfocal
  real,    save :: sim_zc
  real,    save :: sim_zw

  integer, parameter :: sim_refinementLevelMax  =  6
  real,    parameter :: sim_keV2Kelvin          =  11604505.

  integer, save, allocatable :: sim_radialBox (:)
  integer, save, allocatable :: sim_xvalueBox (:)
  real,    save, allocatable :: sim_boxRadius (:)
  real,    save, allocatable :: sim_boxXvalue (:)

end module Simulation_data
