!!****if* source/Simulation/SimulationMain/Cool_Test/Simulation_data
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
!!  Store the simulation data for Globular cluster problem
!!
!! ARGUMENTS
!!
!!
!! PARAMETERS
!!
!!
!!   
!!
!!***

module Simulation_data

  implicit none
#include "constants.h"

  !! These are the physical constants we use
  real, save :: sim_gasConst
  real, save :: sim_gamma

 real, save    :: sim_xMin, sim_xMax, sim_yMin, sim_yMax, sim_zMin, sim_zMax
 real, save    :: sim_c_den, sim_c_temp
 real, save    :: sim_constrast

!!Fluid

    real, save :: sim_xH,sim_xHP,sim_xHM,sim_xD,sim_xDP,sim_xDM
    real, save :: sim_xHE,sim_xHEP,sim_xHEPP,sim_xH2,sim_xH2P
    real, save :: sim_xHD,sim_xHDP,sim_xELEC, sim_xD2, sim_xD2P

!!CHEMCOOL
    real, save :: sim_pchem_time
    real, save :: sim_cool_time

!!METALS

    real, save :: sim_meta
end module Simulation_data
