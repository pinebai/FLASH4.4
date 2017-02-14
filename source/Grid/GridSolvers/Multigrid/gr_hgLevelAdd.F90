!!****if* source/Grid/GridSolvers/Multigrid/gr_hgLevelAdd
!!
!! NAME
!!  gr_hgLevelAdd
!!
!! SYNOPSIS
!!
!!  gr_hgLevelAdd(integer, intent(in) :: level,
!!                integer, intent(in) :: ivar1,
!!                integer, intent(in) :: ivar2,
!!                integer, intent(in) :: LeafFlag)
!!
!! DESCRIPTION
!!
!!  Add the contents of the mesh variable ivar1 to ivar2, placing the result in
!!  ivar1.
!!
!! ARGUMENTS
!!
!!  level    - The refinement of the leaf blocks to add
!!  ivar1    - the first variable index
!!  ivar2    - the second variable index
!!  LeafFlag - can be
!!             MG_NODES_LEAF_ONLY   only add the leaf blocks
!!             MG_NODES_PARENT_ONLY only add the parent blocks
!!             MG_NODES_ALL_NODES   add all blocks
!!
!! RESULT
!!  The variables are added and placed in ivar1
!!
!!***

!!REORDER(5): unk

subroutine gr_hgLevelAdd(level, ivar1, ivar2, LeafFlag)
  use tree, ONLY : lnblocks,lrefine,nodetype
  use physicaldata, ONLY : unk

#include "constants.h"
#include "Multigrid.h"

  implicit none
  
  integer, intent(in) :: ivar1, ivar2, level, LeafFlag

  integer             :: b
  logical             :: AddThisBlock

!==============================================================================


  do b = 1, lnblocks
     
     AddThisBlock = (lrefine(b) == level)
     if (LeafFlag == MG_NODES_LEAF_ONLY) then
        AddThisBlock = (AddThisBlock .and. (nodetype(b) == LEAF))
     else if (LeafFlag == MG_NODES_PARENT_ONLY) then
        AddThisBlock = (AddThisBlock .and. (nodetype(b) /= LEAF))
     endif
     
     if (AddThisBlock) then
        
        unk(ivar1,:,:,:,b) = unk(ivar1,:,:,:,b) + unk(ivar2,:,:,:,b)
        
     endif
     
  enddo
  
  !=================================================================

  return
end subroutine gr_hgLevelAdd
