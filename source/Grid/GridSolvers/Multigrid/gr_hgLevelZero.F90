!!****if* source/Grid/GridSolvers/Multigrid/gr_hgLevelZero
!!
!! NAME
!!  gr_hgLevelZero
!!
!! SYNOPSIS
!!
!!  gr_hgLevelZero(integer, intent(in) :: level,
!!                 integer, intent(in) :: ivar,
!!                 integer, intent(in) :: LeafFlag)
!!
!! DESCRIPTION
!!
!!  Zeroes out certain blocks of a given mesh variable
!!
!! ARGUMENTS
!!
!!  level    - The refinement of the leaf blocks to add
!!  ivar     - the variable index
!!  LeafFlag - can be
!!             MG_NODES_LEAF_ONLY    only zero the leaf blocks
!!             MG_NODES_PARENT_ONLY  only zero the parent blocks
!!             MG_NODES_ALL_NODES    zero all blocks
!!
!! RESULT
!!  The variable is zero'd.
!!
!!***

!!REORDER(5): unk

subroutine gr_hgLevelZero(level, ivar, LeafFlag)

  use tree, ONLY : lrefine,lnblocks,nodetype
  use physicaldata, ONLY: unk
  
  implicit none
  
#include "Multigrid.h"
#include "constants.h"
  
  integer, intent(in) :: ivar, level, LeafFlag
  
  integer             :: b
  logical             :: ZeroThisBlock

!==============================================================================

  
  do b = 1, lnblocks
     
     ZeroThisBlock = (lrefine(b) == level)
     if (LeafFlag == MG_NODES_LEAF_ONLY) then
        ZeroThisBlock = (ZeroThisBlock .and. (nodetype(b) == LEAF))
     else if (LeafFlag == MG_NODES_PARENT_ONLY) then
        ZeroThisBlock = (ZeroThisBlock .and. (nodetype(b) /= LEAF))
     endif
     
     if (ZeroThisBlock) then 
        unk(ivar,:,:,:,b) = 0.0
     endif
     
  enddo
!==============================================================================

return
end subroutine gr_hgLevelZero
