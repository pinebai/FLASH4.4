!!****if* source/Simulation/SimulationMain/unitTest/Pipeline/Simulation_data
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
!!  Stores the local data for the pipeline unit test.
!!  
!!***

Module Simulation_data

  implicit none

#include "constants.h"

  character (len = MAX_STRING_LENGTH), save :: sim_baseName

  integer, save :: sim_channelSize
  integer, save :: sim_cubeSide
  integer, save :: sim_depthPascalTriangle
  integer, save :: sim_globalComm
  integer, save :: sim_globalMe
  integer, save :: sim_globalNumProcs
  integer, save :: sim_itemSize
  integer, save :: sim_lowestNumItemsOnProc
  integer, save :: sim_maxItems
  integer, save :: sim_maxItemsPipeline
  integer, save :: sim_numItems
  integer, save :: sim_numSendingChannels
  integer, save :: sim_storedItems

  integer, save :: sim_sendingChannels (1:10)

  real, allocatable :: sim_items (:,:)

end module Simulation_data
