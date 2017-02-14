!!****if* source/Grid/GridSolvers/Multigrid/gr_hgLevelAddScalar
!!
!! NAME
!!  gr_hgLevelAddScalar
!!
!! SYNOPSIS
!!
!!  call gr_hgLevelAddScalar(integer(in) :: level,
!!                           integer(in) :: ivar,
!!                           real(in)    :: scalar,
!!                           integer(in) :: LeafFlag)
!!
!! DESCRIPTION
!!
!!  Adds a scalar to a subset of the blocks
!!
!! ARGUMENTS
!!
!!  level     - The refinement level of the blocks to add to;
!!              use -1 to not filter by refinement level.
!!  ivar      - the first variable index
!!  scalar    - the scalar to add
!!  LeafFlag  - determines which blocks are modified.
!!              MG_NODES_LEAF_ONLY  only add to the leaf blocks
!!              MG_NODES_PARENT_ONLY  only add to the non-leaf  blocks
!!              MG_NODES_ALL_NODES  add to all blocks
!!
!! RESULT
!!  The variable is uniformly increased by scalar.
!!
!!***

!!REORDER(5): unk

subroutine gr_hgLevelAddScalar(level, ivar, scalar, LeafFlag)

!================================================================

  use tree, ONLY : lrefine,lnblocks,nodetype
  use physicaldata, ONLY: unk

  implicit none

#include "Multigrid.h"
#include "constants.h"

  integer, intent(in) :: ivar, level, LeafFlag
  real, intent(in)    :: scalar

  integer             :: b
  logical             :: AddThisBlock

!==============================================================================

  
  do b = 1, lnblocks
     
     AddThisBlock = (lrefine(b) == level .OR. level == -1)
     if (LeafFlag == MG_NODES_LEAF_ONLY) then
        AddThisBlock = (AddThisBlock .and. (nodetype(b) == LEAF))
     else if (LeafFlag == MG_NODES_PARENT_ONLY) then
        AddThisBlock = (AddThisBlock .and. (nodetype(b) /= LEAF))
     endif
     
     if (AddThisBlock) then
        unk(ivar,:,:,:,b) = unk(ivar,:,:,:,b) + scalar
     
     endif
     
  enddo

  !========================================================================

  return
end subroutine gr_hgLevelAddScalar
