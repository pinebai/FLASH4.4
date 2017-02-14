!!****if* source/Simulation/SimulationMain/unitTest/Grid/GCell/Simulation_data
!!
!! NAME
!!
!!  Simulation_data
!!
!! SYNOPSIS
!!  module Simulation_data()
!!
!! DESCRIPTION
!!
!!  Store the simulation data for unitTesting of guardcell fill
!!   
!! ARGUMENTS
!!
!! PARAMETERS   
!!      Described in the Config file
!!
!! NOTES
!!
!!  No arguments.  All data passed by "use Simulation_data"
!!
!!***

module Simulation_data

  implicit none
#include "constants.h"
#include "Flash.h"
  
  real, save :: sim_smlrho
  integer, save :: sim_subSample


  integer, save :: sim_meshMe
end module Simulation_data
