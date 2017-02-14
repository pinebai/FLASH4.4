!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhFinalizeBlock
!!
!! NAME
!!
!!  gr_bhFinalizeBlock
!!
!!
!! SYNOPSIS
!!
!!  call gr_bhFinalizeBlock(
!!                          integer(in)    :: blockno,
!!                          integer(in)    :: blkLimits(2,MDIM),
!!                          real,pointer   :: solnData(:,:,:,:)
!!        )
!!
!! DESCRIPTION
!!
!!  Called during the Tree Walk after the tree is traversed for all
!!  cells of the block.
!! ARGUMENTS
!!
!!  blockno     : number of block for which the tree walk was just finished
!!  blkLimits   : limits of indeces in the block
!!  solnData    : solution data from the grid
!!
!!
!!***

subroutine gr_bhFinalizeBlock(blockno, blkLimits, solnData)
  use Gravity_interface, ONLY : Gravity_bhFinalizeBlock
  use TreeRay_interface, ONLY : TreeRay_bhFinalizeBlock
  implicit none
#include "constants.h"
  integer, intent(IN) :: blockno
  integer, dimension(2,MDIM), intent(IN)   :: blkLimits
  real, DIMENSION(:,:,:,:), POINTER :: solnData

  ! call _bhFinalizeBlock subroutines of physical modules
  call Gravity_bhFinalizeBlock(blockno, blkLimits, solnData)
  call TreeRay_bhFinalizeBlock(blockno, blkLimits, solnData)
  
  return
end subroutine gr_bhFinalizeBlock


