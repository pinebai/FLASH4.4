!!****if* source/Simulation/SimulationMain/unitTest/RungeKutta/2Dellipse/Simulation_data
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
!!  Stores the local data for the Runge Kutta 2D ellipse unit test.
!!  
!!***

Module Simulation_data

  implicit none

#include "constants.h"
#include "Flash.h"

  character (len = MAX_STRING_LENGTH), save :: sim_RungeKuttaMethod

  integer, save :: sim_numberOfEllipses

  real,    save :: sim_deg2rad
  real,    save :: sim_ellipseAspectRatio
  real,    save :: sim_ellipseCenterX
  real,    save :: sim_ellipseCenterY
  real,    save :: sim_ellipseMajorSemiAxis
  real,    save :: sim_ellipseMinorSemiAxis
  real,    save :: sim_ellipseRotationAngle
  real,    save :: sim_errorFraction
  real,    save :: sim_k
  real,    save :: sim_stepSize
  real,    save :: sim_x0
  real,    save :: sim_y0

end module Simulation_data
