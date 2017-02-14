!!****if* source/Simulation/SimulationMain/unitTest/Pipeline/Cubic/sim_calculateMaxItems
!!
!! NAME 
!!
!!  sim_calculateMaxItems
!!
!! SYNOPSIS
!!
!!  sim_calculateMaxItems ()
!!
!! DESCRIPTION
!!
!!  This function calculates the maximum number of items to be considered.
!!
!! ARGUMENTS
!!
!!  none
!!
!!***

integer function sim_calculateMaxItems ()

  use Simulation_Data,  ONLY : sim_lowestNumItemsOnProc

  implicit none

  sim_calculateMaxItems = sim_lowestNumItemsOnProc

  return
end function sim_calculateMaxItems
