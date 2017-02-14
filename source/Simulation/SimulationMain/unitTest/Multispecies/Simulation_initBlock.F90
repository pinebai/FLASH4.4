!!****if* source/Simulation/SimulationMain/unitTest/Multispecies/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(in)::blockID) 
!!
!!
!!
!! DESCRIPTION
!!
!!  Dummy routine, not needed by Multispecies unit test
!! 
!! ARGUMENTS
!!
!!  blockID:         the number of the block to update
!!  
!!
!!
!!***

subroutine Simulation_initBlock(blockID)

  implicit none
  
  integer, intent(in) :: blockID
  
  
  return
end subroutine Simulation_initBlock
