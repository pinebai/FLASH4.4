!!****if* source/Simulation/SimulationMain/unitTest/Eos/Simulation_data
!!
!! NAME
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  use Simulation_data
!!   
!!
!!***

module Simulation_data
  implicit none
  real, save :: sim_xmin, sim_xmax, sim_ymin, sim_ymax, sim_zmin, sim_zmax

  real, save :: sim_smallx

  integer, save ::  sim_initialMass

  real, save :: sim_densMin, sim_tempMin, sim_xnMin, sim_presMin
  real, save :: sim_densMax, sim_tempMax, sim_xnMax, sim_presMax

integer, save :: sim_meshMe
end module Simulation_data
