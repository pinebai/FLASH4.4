!!****if* source/Simulation/SimulationMain/Plasma/Simulation_data
!!
!! NAME
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  use Simulation_data 
!!
!!  DESCRIPTION
!!
!!  Stores the local data for Simulation setup: Plasma
!!  
!!  
!!
!!***

module Simulation_data
  implicit none
#include "constants.h"
  
  real, save    :: sim_pmass_1
  real, save    :: sim_pvelx_1, sim_pvely_1, sim_pvelz_1

  real, save    :: sim_pmass_2

  real, save    :: sim_bx
  real, save    :: sim_by
  real, save    :: sim_bz

  real, save :: sim_bsw(MDIM), sim_esw(MDIM), sim_usw(MDIM)

end module Simulation_data
