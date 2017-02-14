!!****if* source/Simulation/SimulationMain/StirTurb/Simulation_data
!!
!! NAME
!!
!!  Simulation_data
!!
!! SYNOPSIS
!!   
!!   use Simulation_data
!!
!! DESCRIPTION
!!
!!  Stores the local data for Simulation setup: StirTurb
!!
!!
!!***


module Simulation_data

  implicit none
#include "Eos.h"

  integer, save :: sim_meshMe

  !! *** Runtime Parameters *** !!
  real, save :: sim_smallX
  real, save :: sim_rhoAmbient, sim_cAmbient, sim_mach, sim_gamma

  !! *** EOS Parameters *** !!
  real, save, dimension(EOS_NUM) :: sim_eosArr
  integer, save :: sim_vecLen, sim_mode

end module Simulation_data


