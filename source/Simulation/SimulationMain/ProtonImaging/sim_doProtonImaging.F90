!!****if* source/Simulation/SimulationMain/ProtonImaging/sim_doProtonImaging
!!
!!  NAME 
!!
!!   sim_doProtonImaging
!!
!!  SYNOPSIS
!!
!!   sim_doProtonImaging ()
!!
!!  DESCRIPTION
!!
!!   This routine calls the proton imaging main routine.
!!
!! ARGUMENTS
!!
!!***

subroutine sim_doProtonImaging ()

  use Grid_interface,           ONLY : Grid_getListOfBlocks
  use ProtonImaging_interface,  ONLY : ProtonImaging
  use Driver_data,              ONLY : dr_simTime

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
  call ProtonImaging (blockCount, blockList, dr_simTime)
!
!
!     ...Ready!
!
!
  return
end subroutine sim_doProtonImaging
