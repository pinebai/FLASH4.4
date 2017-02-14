!!****if* source/Simulation/SimulationMain/unitTest/PFFT_BlktriFD/Simulation_data
!!
!! NAME
!!
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  Simulation_data()
!!
!! DESCRIPTION
!!
!!  Stores the local data for Simulation setup: INS-iso-turb
!!
!!***

module Simulation_data

  implicit none

#include "constants.h"

  !! *** Runtime Parameters *** !!
  real, save    :: sim_xMin, sim_xMax, sim_yMin, sim_yMax, sim_zMin, sim_zMax
  logical, save :: sim_gCell
  real, save :: pfb_waven_x
  real, save :: pfb_waven_y
  real, save :: pfb_waven_z
  real, save :: pfb_alpha_x

  integer, save :: sim_meshMe
end module Simulation_data
