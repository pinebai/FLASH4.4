!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhDestroyTree
!!
!! NAME
!!
!!  gr_bhDestroyTree
!!
!!
!! SYNOPSIS
!!
!!   call gr_bhDestroyTree()
!!
!! DESCRIPTION
!!
!!   Deallocates all block-trees.
!!
!! ARGUMENTS
!!
!!
!!***



subroutine gr_bhDestroyTree()

  use Grid_interface, ONLY : Grid_getLocalNumBlks, &
      Grid_getListOfBlocks, Grid_updateRefinement
  use gr_bhData, ONLY : gr_bhTreeNumProcs, gr_bhTreeNodetype, &
    gr_bhTreeArray
  use Timers_interface, ONLY : Timers_start, Timers_stop
      
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  integer :: localNumBlocks, i
  integer :: blockID, blockCount
  integer :: blockList(MAXBLOCKS)

  call Timers_start("destroy_tree")

  do i = 0, gr_bhTreeNumProcs-1
    do blockID = 1, MAXBLOCKS !gr_bhTreeNumBlocks(i)
      if (gr_bhTreeNodetype(blockID, i) .ne. 1) cycle
      deallocate(gr_bhTreeArray(i, blockID)%p)
    enddo
  enddo

  call Timers_stop("destroy_tree")
  return
end subroutine gr_bhDestroyTree


