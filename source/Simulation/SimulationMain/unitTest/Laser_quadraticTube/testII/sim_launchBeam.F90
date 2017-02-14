!!****if* source/Simulation/SimulationMain/unitTest/Laser_quadraticTube/testII/sim_launchBeam
!!
!!  NAME 
!!
!!   sim_launchBeam
!!
!!  SYNOPSIS
!!
!!   sim_launchBeam ()
!!
!!  DESCRIPTION
!!
!!   This routine launches 1 circular beam of rays onto the xz-plane of a quadratic
!!   potential tube and records their (x,z) exit points and their exit power. The analysis of
!!   these results is done in a separate routine.
!!
!! ARGUMENTS
!!
!!***

subroutine sim_launchBeam ()

  use Grid_interface,              ONLY : Grid_getListOfBlocks

  use EnergyDeposition_interface,  ONLY : EnergyDeposition

  implicit none

#include "Flash.h"
#include "constants.h"

  integer  :: blockCount
  integer  :: pass

  real     :: dt
  real     :: time

  integer  :: blockList (1:MAXBLOCKS)
!
!
!     ...Get the list of blocks on current processor.
!
!
  call Grid_getListOfBlocks (LEAF,   blockList, blockCount)
!
!
!     ...Launch all rays in beam by calling the energy deposition routine.
!
!
  dt   = 1.e-8      ! this ensures correct initial ray power of 1 erg/sec
  pass = 1          ! mimics unsplit energy deposition
  time = 0.0        ! overall simulation time

  call EnergyDeposition (blockCount, blockList, dt, time, pass)
!
!
!     ...Ready!
!
!
  return
end subroutine sim_launchBeam
