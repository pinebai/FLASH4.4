!!****if* source/Simulation/SimulationMain/unitTest/RungeKutta/3Dcircle/Simulation_data
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
!!  Stores the local data for the Runge Kutta unit test.
!!  
!!***

Module Simulation_data

  implicit none

#include "constants.h"
#include "Flash.h"

  character (len = MAX_STRING_LENGTH), save :: sim_RungeKuttaMethod

  integer, save :: sim_numberOfCircles
  integer, save :: sim_numberOfRungeKuttaSteps

  real,    save :: sim_errorFraction
  real,    save :: sim_radius
  real,    save :: sim_speed
  real,    save :: sim_stepSize

  real,    save :: sim_rx0
  real,    save :: sim_ry0
  real,    save :: sim_rz0

  real,    save :: sim_vx0
  real,    save :: sim_vy0
  real,    save :: sim_vz0

end module Simulation_data
