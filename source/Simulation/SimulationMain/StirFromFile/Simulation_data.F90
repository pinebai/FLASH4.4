!!****if* source/Simulation/SimulationMain/StirFromFile/Simulation_data
!!
!! NAME
!!  Simulation_data
!!
!! SYNOPSIS
!!   use Simulation_data
!!
!! DESCRIPTION
!!  Stores the local data for Simulation setup: StirFromFile
!!
!!***


module Simulation_data

  implicit none

  !! *** Runtime Parameters *** !!
  real, save    :: sim_rhoAmbient, sim_cAmbient, sim_gamma, sim_MagField_z
  logical, save :: sim_magnetic

end module Simulation_data


