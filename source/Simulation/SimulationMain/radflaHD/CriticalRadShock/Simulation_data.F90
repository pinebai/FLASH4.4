!!****if* source/Simulation/SimulationMain/radflaHD/CriticalRadShock/Simulation_data
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
!!  Stores the simulation data
!!
!! 
!!***

module Simulation_data

  implicit none
  
#include "constants.h"

  !! *** Runtime Parameters *** !!
  real, save :: sim_rho
  real, save :: sim_tgas
  real, save :: sim_trad
  real, save :: sim_velx
  real, save :: sim_xmin, sim_xmax
  real, save :: sim_ymin, sim_ymax
  real, save :: sim_zmin, sim_zmax
  real, save :: sim_gamma

  logical, save :: sim_gCell = .true., sim_no_gCell = .false.
 
end module Simulation_data


