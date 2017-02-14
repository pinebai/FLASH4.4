!!****if* source/Simulation/SimulationMain/RHD_Riemann2D/Simulation_data
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
!!  Store the simulation data for the Sod problem
!!
!! ARGUMENTS
!!
!!  None
!!
!!***

module Simulation_data

  implicit none

  real, save :: sim_gamma, sim_smallP, sim_smallX
  logical, save :: sim_gCell

  integer, save :: sim_meshMe
end module Simulation_data


