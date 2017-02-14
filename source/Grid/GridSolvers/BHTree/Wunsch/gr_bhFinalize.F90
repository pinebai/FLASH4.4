!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhFinalize
!!
!! NAME
!!
!!  gr_bhFinalize
!!
!! 
!! SYNOPSIS
!!
!!  call gr_bhFinalize()
!!
!!
!! DESCRIPTION
!!
!!  Finalizes the tree Poisson solver.  Deallocates all arrays.
!!
!!***

subroutine gr_bhFinalize()

  use gr_bhData, ONLY: gr_bhTreeArray, &
       gr_bhLocSentTreeLevels,gr_bhLocRecvTreeLevels,gr_bhTreeCellSize, &
       gr_bhTreeNodeSize,gr_bhTreeNodeSize2,gr_bhTreeBCen,gr_bhLocCoords, &
       gr_bhTreeParentTree,gr_bhLocParentTree, &
       gr_bhTreeNodetype,gr_bhTreeLrefine,gr_bhTreeChild,gr_bhTreeLnblocks, &
       gr_bhTreeBlocklist, gr_bhTreeFirstLevBlocks, gr_bhBlockTreePos, &
       gr_bhPriorityQueue

  implicit none

  deallocate(gr_bhTreeArray)
  deallocate(gr_bhLocSentTreeLevels)
  deallocate(gr_bhLocRecvTreeLevels)
  deallocate(gr_bhTreeCellSize)
  deallocate(gr_bhTreeNodeSize)
  deallocate(gr_bhTreeNodeSize2)
  deallocate(gr_bhTreeBCen)
  deallocate(gr_bhLocCoords)
  deallocate(gr_bhBlockTreePos)

  deallocate(gr_bhTreeParentTree)
  deallocate(gr_bhLocParentTree)

  deallocate(gr_bhTreeNodetype)
  deallocate(gr_bhTreeLrefine)
  deallocate(gr_bhTreeChild)
  deallocate(gr_bhTreeLnblocks)
  deallocate(gr_bhTreeBlocklist)
  deallocate(gr_bhTreeFirstLevBlocks)

  deallocate(gr_bhPriorityQueue)

  return
end subroutine gr_bhFinalize
