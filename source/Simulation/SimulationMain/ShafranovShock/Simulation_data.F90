!!****if* source/Simulation/SimulationMain/ShafranovShock/Simulation_data
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
!!  Store the simulation data for the Shafranov problem
!!
!! ARGUMENTS
!!
!!
!! PARAMETERS
!!
!!  sim_DataPoint : The number of data points in init file.   
!!  sim_InitData  : The problem is initialized using data from this file.
!!
!!   
!!
!!***

module Simulation_data

  implicit none

#include "constants.h"

  !! *** Runtime Parameters *** !!
  integer,                          save :: sim_DataPoints
  character(len=MAX_STRING_LENGTH), save :: sim_InitData
  real                            , save :: sim_ShockSpeed
  real,                             save :: sim_smallX
  real,                             save :: sim_maxTol, sim_abar, sim_zbar

  !! *** Variables pertaining to Simulation Setup 'Shafranov' *** !!
  real, allocatable,dimension(:),   save :: sim_x, sim_velx, sim_tele, sim_rho, &
                                            sim_tion, sim_pele, sim_pion,       &
                                            sim_sele

  integer, save :: sim_meshMe
end module Simulation_data


