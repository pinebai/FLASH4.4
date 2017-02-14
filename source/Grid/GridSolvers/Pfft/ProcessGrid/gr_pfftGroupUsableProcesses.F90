!!****if* source/Grid/GridSolvers/Pfft/ProcessGrid/gr_pfftGroupUsableProcesses
!!
!! NAME
!!
!!  gr_pfftGroupUsableProcesses
!!
!! SYNOPSIS
!!  
!!  gr_pfftGroupUsableProcesses(integer(IN)  :: myPE,
!!                              integer(IN)  :: globalProcs
!!                              integer(IN)  :: originalComm,
!!                              integer(OUT)  :: newComm)
!!
!! DESCRIPTION
!!
!!  We specify a communicator which contains each processor 
!!  that we will use in the PFFT framework.  If we have more
!!  blocks than processors then we can safely use every single 
!!  processor, otherwise we must create a subset of processors.
!!
!! ARGUMENTS
!!
!!  myPE - my rank in FLASH grid.
!!  globalProcs - total number of processors available.
!!  originalComm - the original communicator.
!!  newComm - the new communicator.
!!
!!***
subroutine gr_pfftGroupUsableProcesses(myPE, globalProcs, originalComm, newComm)

#include "Flash.h"
#include "constants.h"

#ifdef FLASH_GRID_PARAMESH
  use tree, ONLY : lnBlocks, lrefine
  use Grid_data, ONLY : gr_oneRefLev, gr_nblockX, gr_nblockY, gr_nblockZ
#endif
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none
  integer, intent(IN) :: myPE, globalProcs, originalComm
  integer, intent(OUT) :: newComm
  integer, dimension(1:MDIM) :: initialBlocks
  integer :: totalBlocks, numBlocksAlongDimension, communicatorIndicator, &
       i, ierr, refinementLevel

#include "Flash_mpi.h"

#if defined(FLASH_GRID_UG)

  !By definition we have one block per process.
  newComm = originalComm
  return

#elif defined(FLASH_GRID_PARAMESH)

  refinementLevel = gr_oneRefLev
  numBlocksAlongDimension = 2 ** (refinementLevel-1)
  initialBlocks(IAXIS) = gr_nblockX
  initialBlocks(JAXIS) = gr_nblockY
  initialBlocks(KAXIS) = gr_nblockZ

  !Total block existing in FLASH grid at this refinement level.
  totalBlocks = product(initialBlocks(1:NDIM) * numBlocksAlongDimension)


#if defined (DEBUG_PFFT)
  if (myPE == 0) then
     print *, "Processor:", myPE, "calculates", totalBlocks, &
          "usable blocks at refinement level:", refinementLevel
  end if
#endif


  !Check if there are more blocks at this refinement level than the
  !global number of processors.
  if (totalBlocks < globalProcs) then
     if (myPE == 0) then
        print *, "[Warning]: More processors than blocks, so using process subset"
     end if

     communicatorIndicator = MPI_UNDEFINED   !Important.

     !Check each block on this processor, and specify that this processor
     !will be part of the communicator if it has at least one block at this
     !refinement level.
     do i = 1, lnBlocks
        if (lrefine(i) == refinementLevel) then
           communicatorIndicator = 1
           exit
        end if
     end do

     call MPI_Comm_split(originalComm, communicatorIndicator, myPE, &
          newComm, ierr)

  else

     newComm = originalComm  !We are OK to use all available processors.

  end if

#endif

end subroutine gr_pfftGroupUsableProcesses
