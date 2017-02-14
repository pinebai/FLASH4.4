!!****if* source/Simulation/SimulationMain/magnetoHD/BrioWu/Simulation_data
!!
!! NAME
!!
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  use Simulation_data
!!
!! DESCRIPTION
!!
!!  Stores the local data for Simulation setup: BrioWu
!!
!!***

module Simulation_data

  implicit none

#include "constants.h"

  !! *** Runtime Parameters *** !!    
  real, save   :: sim_posn, sim_gamma
  real, save   :: sim_uRight, sim_uLeft, sim_vRight, sim_vLeft, sim_wLeft, sim_wRight
  real, save   :: sim_rhoRight, sim_rhoLeft, sim_pRight, sim_pLeft
  real, save   :: sim_byLeft, sim_byRight, sim_bzLeft, sim_bzRight, sim_bNormal
  real, save   :: sim_xmin,sim_xmax,sim_ymin,sim_ymax

  !! Simulation variables
  real, save    :: sim_xangle, sim_yangle, sim_xcos, sim_ycos, sim_zcos
  real, save    :: sim_smallx, sim_smallP

  logical, save :: sim_gCell, sim_killdivb

  integer, save :: sim_meshMe
end module Simulation_data
