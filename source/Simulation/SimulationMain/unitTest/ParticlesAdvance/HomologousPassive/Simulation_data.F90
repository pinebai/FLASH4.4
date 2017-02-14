!!****if* source/Simulation/SimulationMain/unitTest/ParticlesAdvance/HomologousPassive/Simulation_data
!!
!! NAME
!!
!!  Simulation_data
!!
!! SYNOPSIS
!!  use Simulation_data
!!
!! DESCRIPTION
!!
!!  Store the simulation data for unitTesting of Particles
!!   
!! ARGUMENTS
!!
!! PARAMETERS   
!!      Described in the Config file
!!
!!
!!***

module Simulation_data

  implicit none
#include "constants.h"
#include "Flash.h"
  
  real, save    :: sim_rho_amb, sim_p_amb, sim_vx_amb, sim_vx_multiplier
  real, save    :: sim_seed, sim_vx_pert, sim_vy_pert, sim_vz_pert
  real, save    :: sim_a0, sim_a1

  logical,save :: sim_analyticParticlePositions, sim_fakeMapMeshToParticles

  integer,save :: sim_schemeOrder
  real,save :: sim_maxTolCoeff0,sim_maxTolCoeff1,sim_maxTolCoeff2,sim_maxTolCoeff3

integer, save :: sim_meshMe
end module Simulation_data
