!!****if* source/Simulation/SimulationMain/SodSpherical/Simulation_data
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
!!  Store the simulation data for the SodSpherical problem
!!
!! ARGUMENTS
!!
!!
!! PARAMETERS
!!
!!  sim_rhoLeft    Density in the left part of the grid
!!  sim_rhoRight   Density in the right part of the grid
!!  sim_pLeft      Pressure  in the left part of the grid
!!  sim_pRight     Pressure  in the righ part of the grid
!!  sim_uLeft      fluid velocity in the left part of the grid
!!  sim_uRight     fluid velocity in the right part of the grid
!!  sim_xangle     Angle made by diaphragm normal w/x-axis (deg)
!!  sim_yangle     Angle made by diaphragm normal w/y-axis (deg)
!!  sim_shockpos   Point of intersection between the shock plane and the x-axis
!!
!!
!!   
!!
!!***

module Simulation_data

  implicit none

  !! *** Runtime Parameters *** !!

  real, save :: sim_rhoLeft, sim_rhoRight, sim_pLeft, sim_pRight
  integer, save :: sim_idir
  real, save :: sim_shockpos
  real, save :: sim_gamma, sim_smallX
  real, save :: sim_hydroCfl = 0.0

integer, save :: sim_meshMe
end module Simulation_data


