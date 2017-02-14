!!****if* source/Simulation/SimulationMain/Pancake/Simulation_data
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
!!  Stores the local data for Simulation setup:  Pancake
!!  
!!
!!***

module Simulation_data

  implicit none

#include "constants.h"

  integer, save :: sim_meshMe
  
  !! *** Runtime Parameters *** !!
  
  real, save    :: sim_xMin, sim_xMax, sim_yMin, sim_yMax, sim_zMin, sim_zMax
  real, save    :: sim_hubble, sim_lambda, sim_omegaMatter, sim_omegaBaryon
  real, save    :: sim_Tfiducial, sim_zfiducial, sim_zcaustic, sim_wavelength
  real, save    :: sim_xangle, sim_yangle
  real, save    :: sim_smallE, sim_smallRho, sim_smallP, sim_gamma
  real, save    :: sim_Newton, sim_pi, sim_gascon

  logical, save :: sim_useParticles

  integer, save :: sim_numPEs, sim_nxp, sim_nyp, sim_nzp
  real, save    :: sim_mr

  !! *** Auxiliary Variables - Introduced in Simulation_init *** !!

  real, save    :: sim_zinitial, sim_critden, sim_meanden, sim_bmeanden
  real, save    :: sim_dmeanden, sim_kmag, sim_xcos, sim_ycos, sim_zcos

  real, save    :: sim_temp1, sim_temp2, sim_temp3

end module Simulation_data
