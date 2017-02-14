!!****if* source/Simulation/SimulationMain/unitTest/PhysConst/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(blockID)
!!
!! DESCRIPTION
!!
!!  The unitTest for Physical Constants does not use any physical domain.
!!    Therefore, this routine is included only to maintain the API for
!!    Simulation units.
!!
!!***

subroutine Simulation_initBlock(blockID)
#include "constants.h"

  implicit none
  
  integer, intent(in) :: blockID
  
  
  return
end subroutine Simulation_initBlock
