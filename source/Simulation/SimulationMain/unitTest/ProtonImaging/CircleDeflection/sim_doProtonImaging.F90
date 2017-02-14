!!****if* source/Simulation/SimulationMain/unitTest/ProtonImaging/CircleDeflection/sim_doProtonImaging
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
!!   This routine performs the proton imaging.
!!
!! ARGUMENTS
!!
!!***

subroutine sim_doProtonImaging ()

  use Grid_interface,           ONLY : Grid_getListOfBlocks

  use ProtonImaging_interface,  ONLY : ProtonImaging

  implicit none

#include "Flash.h"
#include "constants.h"

  integer  :: blockCount

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
!     ...Call the main proton imaging routine.
!
!
  time = 0.0        ! overall simulation time

  call ProtonImaging (blockCount, blockList, time)
!
!
!     ...Ready!
!
!
  return
end subroutine sim_doProtonImaging
