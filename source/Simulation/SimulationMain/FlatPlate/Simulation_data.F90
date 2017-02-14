!!****if* source/Simulation/SimulationMain/FlatPlate/Simulation_data
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
!!  Store the simulation data for the wind tunnel problem with a step
!!
!! ARGUMENTS
!!
!!
!! PARAMETERS
!!
!!  sim_pAmbient    Initial ambient pressure
!!  sim_rhoAmbient  Initial ambient density
!!  sim_windVel     Inflow velocity (parallel to x-axis)
!!  gamma           the Gamma EOS thing
!!  smallp          minimum for pressure
!!  smallx          minimum for abundance
!!
!!
!!   
!!
!!***

module Simulation_data

  implicit none

  !! *** Runtime Parameters *** !!

  real, save :: sim_pAmbient, sim_rhoAmbient, sim_windVelx, sim_windVely, sim_windVelz, sim_rhoBulk
  real, save :: sim_gamma, sim_smallP, sim_smallX
  real, save :: sim_radius, sim_xCtr, sim_yCtr, sim_zCtr
  real, save :: sim_xCos, sim_yCos, sim_zCos,sim_xAngle, sim_yAngle, sim_zAngle
  real, save :: sim_Mach
  logical, save :: sim_gCell
  integer, save :: sim_meshMe, sim_number


end module Simulation_data


