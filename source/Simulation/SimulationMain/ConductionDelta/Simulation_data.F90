!!****if* source/Simulation/SimulationMain/ConductionDelta/Simulation_data
!!
!! NAME
!!  Simulation_data
!!
!! SYNOPSIS
!!  Use Simulation_data
!!
!! DESCRIPTION
!!
!!  Store the simulation data
!!
!! 
!!***

module Simulation_data

  implicit none

  !! *** Runtime Parameters *** !!

  real, save :: sim_xCenter, sim_yCenter, sim_zCenter
  real, save :: sim_gamma, sim_smallP, sim_smallX
  integer, save ::  sim_orientation
  real, save :: sim_rhoInit, sim_toffset, sim_Q, sim_tempBackground

  logical, save :: sim_gCell

  real, save :: sim_iniCondTemperatureExponent
  real, save :: sim_condTemperatureExponent, sim_alpha, sim_xi0, sim_xfInitial

integer, save :: sim_meshMe
end module Simulation_data


