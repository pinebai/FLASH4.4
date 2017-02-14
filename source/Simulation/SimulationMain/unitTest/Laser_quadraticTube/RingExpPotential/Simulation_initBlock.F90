!!****if* source/Simulation/SimulationMain/unitTest/Laser_quadraticTube/RingExpPotential/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock (integer (in) :: blockID)
!!
!! DESCRIPTION
!!
!!  Initializes the data for a specified block needed to do the laser quadratic tube ring problem.
!!
!! ARGUMENTS
!!
!!  blockID : The block ID number to be initialized
!!
!!***

subroutine Simulation_initBlock (blockID)

  use Driver_interface,  ONLY : Driver_abortFlash

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent (in) :: blockID
!
!
!    ...Call the 3D routine.
!
!
  call sim_initBlock3DRec (blockID)
!
!
!    ...Ready!
!
!
  return
end subroutine Simulation_initBlock
