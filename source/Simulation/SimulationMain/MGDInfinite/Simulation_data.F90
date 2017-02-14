!!****if* source/Simulation/SimulationMain/MGDInfinite/Simulation_data
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

#include "Flash.h"

  implicit none

  integer :: sim_fileUnit
  integer :: sim_globalMe

  !! *** Runtime Parameters *** !!
  real, save :: sim_rho
  real, save :: sim_tion
  real, save :: sim_tele
  real, save :: sim_trad

  real, save :: sim_massfracs(NSPECIES)

  real, save :: sim_smallX

end module Simulation_data


