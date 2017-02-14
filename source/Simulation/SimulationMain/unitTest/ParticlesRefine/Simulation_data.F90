!!****if* source/Simulation/SimulationMain/unitTest/ParticlesRefine/Simulation_data
!!
!! NAME
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  use Simulation_data
!!
!! DESCRIPTION
!!
!!  Store the simulation data for the poistest problem
!!
!! ARGUMENTS
!!
!!
!! PARAMETERS
!!
!!
!!   
!!
!!***

module Simulation_data

  implicit none

  !! *** Runtime Parameters *** !!

  real, save :: sim_smlRho
  real, save :: sim_ptMass, sim_densityThreshold
  integer, save :: sim_meshMe, sim_minBlks
  logical, save :: sim_print=.false.
end module Simulation_data
