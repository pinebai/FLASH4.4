!!****if* source/Simulation/SimulationMain/unitTest/Laser_quadraticTube/RingExpPotential/Simulation_data
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

  integer, save :: sim_alpha
  integer, save :: sim_globalComm
  integer, save :: sim_globalMe
  integer, save :: sim_globalNumProcs
  integer, save :: sim_numRaysAnalyzed
  integer, save :: sim_numRaysLaunched
  integer, save :: sim_refinementLevel
  integer, save :: sim_totalNumberOfBoxes

  real,    save :: sim_beamFrequency
  real,    save :: sim_beamWavelength
  real,    save :: sim_Boltzmann
  real,    save :: sim_cellSizeX
  real,    save :: sim_cellSizeY
  real,    save :: sim_cellSizeZ
  real,    save :: sim_exitPowerFraction
  real,    save :: sim_focalPointRadius
  real,    save :: sim_maxDeltaRsqr
  real,    save :: sim_avgDeltaR
  real,    save :: sim_rmsDeltaR
  real,    save :: sim_sigmaDeltaR
  real,    save :: sim_avgDeltaX
  real,    save :: sim_rmsDeltaX
  real,    save :: sim_sigmaDeltaX

  real,    save :: sim_vy
  real,    save :: sim_nc
  real,    save :: sim_r0
  real,    save :: sim_ne0
  real,    save :: sim_ner0
  real,    save :: sim_Te0
  real,    save :: sim_viB0
  real,    save :: sim_xTubeCenter
  real,    save :: sim_yfocal
  real,    save :: sim_zTubeCenter

  integer, parameter :: sim_refinementLevelMax  =  6
  real,    parameter :: sim_keV2Kelvin          =  11604505.

  real,    save, dimension (1:20) :: sim_FocalTime   = (/    2.00000000000000000   , &
                                                             1.57079632679489662   , &
                                                             1.40218210532545426   , &
                                                             1.31102877714605991   , &
                                                             1.25373062481720727   , &
                                                             1.21432532394379081   , &
                                                             1.18554403060175695   , &
                                                             1.16359257121826938   , &
                                                             1.14629395408666632   , &
                                                             1.13230869752157537   , &
                                                             1.12076700391887598   , &
                                                             1.11107930198320723   , &
                                                             1.10283171867822819   , &
                                                             1.09572512260053927   , &
                                                             1.08953791366386869   , &
                                                             1.08410241908920600   , &
                                                             1.07928942192362020   , &
                                                             1.07499772924602362   , &
                                                             1.07114696271179879   , &
                                                             1.06767246662400211    /)

  real,    save, dimension (1:20) :: sim_IaaIntegral = (/    1.33333333333333333   , &
                                                             7.85398163397448310e-1, &
                                                             5.60872842130181704e-1, &
                                                             4.37009592382019968e-1, &
                                                             3.58208749947773505e-1, &
                                                             3.03581330985947701e-1, &
                                                             2.63454229022612655e-1, &
                                                             2.32718514243653875e-1, &
                                                             2.08417082561212058e-1, &
                                                             1.88718116253595895e-1, &
                                                             1.72425692910596304e-1, &
                                                             1.58725614569029605e-1, &
                                                             1.47044229157097092e-1, &
                                                             1.36965640325067408e-1, &
                                                             1.28180931019278670e-1, &
                                                             1.20455824343245111e-1, &
                                                             1.13609412834065284e-1, &
                                                             1.07499772924602362e-1, &
                                                             1.02013996448742742e-1, &
                                                             9.70611333294547376e-2 /)
  integer, save, allocatable :: sim_powersBox (:)
  integer, save, allocatable :: sim_radialBox (:)
  integer, save, allocatable :: sim_xvalueBox (:)
  real,    save, allocatable :: sim_boxPowers (:)
  real,    save, allocatable :: sim_boxRadius (:)
  real,    save, allocatable :: sim_boxXvalue (:)

end module Simulation_data
