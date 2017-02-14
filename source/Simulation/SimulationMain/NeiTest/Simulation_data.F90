!!****if* source/Simulation/SimulationMain/NeiTest/Simulation_data
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
!!  Stores the local data for Simulation setup: Neitest
!!  
!! PARAMETERS
!!
!!  sim_rhoAmbient     Initial ambient density
!!  sim_velInit        Initial velocity
!!  sim_tAmbient       Initial ambient temperature
!!  sim_tPerturb       Perturbation temperature
!!  sim_radius         Radial position of inner edge of grid (for 1D )
!!
!!
!!***

module Simulation_data

  implicit none
  
  ! save the parameters that describe this initialization
  real, save :: sim_velInit, sim_rhoAmbient, sim_tAmbient, sim_tPerturb
  real, save :: sim_radius
  real, save :: sim_xstep, sim_smallx
  real, save :: vel_init  
  
end module Simulation_data
