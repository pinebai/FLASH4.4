!!****if* source/Grid/GridSolvers/Multigrid/gr_hgLevelMultiplyScalar
!!
!! NAME
!!  gr_hgLevelMultiplyScalar
!!
!! SYNOPSIS
!!
!!  gr_hgLevelMultiplyScalar(integer, intent(in) :: level,
!!                           integer, intent(in) :: ivar,
!!                           real, intent(in)    :: scalar,
!!                           integer, intent(in) :: LeafFlag)
!!
!! DESCRIPTION
!!
!!  Multiply the contents of mesh variable ivar by scalar in-place.
!!
!! ARGUMENTS
!!
!!  level    - The refinement of the leaf blocks to add
!!  ivar     - the first variable index
!!  scalar   - the scalar to multiply by
!!  LeafFlag - can be
!!             MG_NODES_LEAF_ONLY   only add the leaf blocks
!!             MG_NODES_PARENT_ONLY only add the parent blocks
!!             MG_NODES_ALL_NODES   add all blocks
!!
!! RESULT
!!  The variable is multiplied by a factor of ivar
!!
!!***

!!REORDER(5): unk

subroutine gr_hgLevelMultiplyScalar(level, ivar, scalar, LeafFlag)

  use tree, ONLY : lrefine,lnblocks,nodetype
  use physicaldata, ONLY: unk
  
  implicit none

#include "Multigrid.h"  
#include "constants.h"
  
  integer, intent(in) :: ivar, level, LeafFlag
  real, intent(in)    :: scalar
  
  integer             :: b
  logical             :: MultThisBlock

!==============================================================================

  
  do b = 1, lnblocks
     
     MultThisBlock = (lrefine(b) == level)
     if (LeafFlag == MG_NODES_LEAF_ONLY) then
        MultThisBlock = (MultThisBlock .and. (nodetype(b) == LEAF))
     else if (LeafFlag == MG_NODES_PARENT_ONLY) then
        MultThisBlock = (MultThisBlock .and. (nodetype(b) /= LEAF))
     endif
     
     if (MultThisBlock) then 
        unk(ivar,:,:,:,b) = unk(ivar,:,:,:,b) * scalar
     endif
     
  enddo
  
  !========================================================================
  
  return
end subroutine gr_hgLevelMultiplyScalar
