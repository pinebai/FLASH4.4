!!****if* source/Simulation/SimulationMain/Jeans/Simulation_data
!!
!! NAME
!!
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  use Simulation_data 
!!
!!  DESCRIPTION
!!
!!  Stores the local data for Simulation setup: Jeans
!!  
!! PARAMETERS
!!
!!  sim_p0       Initial ambient pressure
!!  sim_rho0     Initial ambient density
!!  sim_lambdaX  Perturbation X-wavelength
!!  sim_lambdaY  Perturbation Y-wavelength
!!  sim_lambdaZ  Perturbation Z-wavelength
!!  sim_A        Perturbation amplitude
!!
!!***

module Simulation_data

  implicit none

#include "constants.h"

  !! *** Runtime Parameters *** !!

  real, save    :: sim_p0, sim_rho0
  real, save    :: sim_gamma, sim_lambdaX, sim_lambdaY, sim_lambdaZ
  real, save    :: sim_A
  real, save    :: sim_smallE
  real, save    :: sim_deltaRef, sim_deltaDeRef, sim_refDensity
  real, save    :: sim_U0 = 0.0 !! thermal energy at time zero

  !! *** Variables pertaining to this Simulation *** !!

  real, save    :: sim_velA, sim_kX, sim_kY, sim_kZ
  real, save    :: sim_kMag, sim_kJ, sim_c0, sim_oscFreq

  !! *** Physical Constants *** !!

  real, save    :: sim_pi, sim_Newton
         
  integer, save :: sim_meshMe
end module Simulation_data
