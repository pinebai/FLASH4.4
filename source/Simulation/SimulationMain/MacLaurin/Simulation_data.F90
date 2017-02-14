!!****if* source/Simulation/SimulationMain/MacLaurin/Simulation_data
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
!!  Stores the local data for Simulation setup:  MacLaurin
!!  
!!
!!***

module Simulation_data

  implicit none

#include "constants.h"
  
  !! *** Runtime Parameters *** !!
  
  real, save    :: sim_eccentricity
  real, save    :: sim_gamma, sim_density, sim_Omega1, sim_a1
  real, save    :: sim_xctr, sim_yctr, sim_zctr
  real, save    :: sim_smallX, sim_smallE, sim_smallRho, sim_smallP
  real, save    :: sim_Newton, sim_pi

  integer, save :: sim_nsubzones, sim_initGeometry
  integer, parameter :: sim_geom2DAxisymmetric=1, sim_geom3DCartesian=2

  !! *** Auxiliary Variables - Introduced in Simulation_init *** !!

  real, save    :: sim_nsubinv, sim_a3, sim_a1inv, sim_a3inv
  real, save    :: sim_Omega2, sim_Pconst


  integer, save :: sim_meshMe
end module Simulation_data
