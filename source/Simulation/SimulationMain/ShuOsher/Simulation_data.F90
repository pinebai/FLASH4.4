!!****if* source/Simulation/SimulationMain/ShuOsher/Simulation_data
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
!!  Store the simulation data for the Shu-Osher problem
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
!!  sim_posnR      Point of intersection between the shock plane and the x-axis
!!
!!
!!***

module Simulation_data

  implicit none

  integer, save :: sim_meshMe

  !! *** Runtime Parameters *** !!

  real, save :: sim_rhoLeft, sim_rhoRight, sim_pLeft, sim_pRight
  real, save :: sim_uLeft, sim_uRight, sim_posn
  real, save :: sim_gamma, sim_smallP, sim_smallX
  real, save :: sim_aRho,  sim_fRho


  !! *** Variables pertaining to Simulation Setup  *** !!
  logical, save :: sim_gCell

end module Simulation_data


