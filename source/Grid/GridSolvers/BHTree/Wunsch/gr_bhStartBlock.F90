!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhStartBlock
!!
!! NAME
!!
!!  gr_bhStartBlock
!!
!!
!! SYNOPSIS
!!
!!  call gr_bhStartBlock(
!!                          integer(in)    :: blockno,
!!                          integer(in)    :: blkLimits(2,MDIM),
!!                          real,pointer   :: solnData(:,:,:,:)
!!        )
!!
!! DESCRIPTION
!!
!!  Called during the Tree Walk before the tree is traversed for a given block.
!!
!! ARGUMENTS
!!
!!  blockno     : number of block for which the tree walk was just finished
!!  blkLimits   : limits of indeces in the block
!!  solnData    : solution data from the grid
!!
!!
!!***

subroutine gr_bhStartBlock(blockno, blkLimits, solnData)
  use Gravity_interface, ONLY : Gravity_bhStartBlock
  use TreeRay_interface, ONLY : TreeRay_bhStartBlock
  implicit none
#include "FortranLangFeatures.fh"
#include "constants.h"
  integer, intent(IN) :: blockno
  integer, dimension(2,MDIM), intent(IN)   :: blkLimits
  real, DIMENSION(:,:,:,:), POINTER_INTENT_IN :: solnData

  ! call _bhStartBlock subroutines of physical modules
  call Gravity_bhStartBlock(blockno, blkLimits, solnData)
  call TreeRay_bhStartBlock(blockno, blkLimits, solnData)
  
  return
end subroutine gr_bhStartBlock


