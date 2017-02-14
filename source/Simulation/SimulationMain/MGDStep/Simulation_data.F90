!!****if* source/Simulation/SimulationMain/MGDStep/Simulation_data
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

#include "constants.h"

  !! *** Runtime Parameters *** !!
  real, save :: sim_rho1
  real, save :: sim_tion1
  real, save :: sim_tele1
  real, save :: sim_trad1

  real, save :: sim_rho2
  real, save :: sim_tion2
  real, save :: sim_tele2
  real, save :: sim_trad2

  !! *** Other Variables *** !!
  character(len=MAX_STRING_LENGTH) :: sim_initGeom
  real, save :: sim_thickness
end module Simulation_data


