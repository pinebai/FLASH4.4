!!****f* source/physics/Gravity/Gravity_bhStartBlock
!!
!! NAME
!!
!!  Gravity_bhStartBlock
!!
!!
!! SYNOPSIS
!!
!!   call Gravity_bhStartBlock(
!!                          integer(in)    :: blockno,
!!                          integer(in)    :: blkLimits(2,MDIM),
!!                          real,pointer   :: solnData(:,:,:,:)
!!        )
!!
!! DESCRIPTION
!!
!!  Called during the Tree Walk before the tree is traversed for a given block.
!!  In this version, calculates the inverted mean error in the acceleration 
!!  for this block.
!!
!! ARGUMENTS
!!
!!  blockno     : number of block for which the tree walk was just finished
!!  blkLimits   : limits of indeces in the block
!!  solnData    : solution data from the grid
!!
!!
!!***

subroutine Gravity_bhStartBlock(blockno, blkLimits, solnData)
  implicit none
#include "constants.h"
#include "Flash.h"
#include "FortranLangFeatures.fh"
  integer, intent(IN) :: blockno
  integer, dimension(2,MDIM), intent(IN)   :: blkLimits
  real, DIMENSION(:,:,:,:), POINTER_INTENT_IN :: solnData

  return
end subroutine Gravity_bhStartBlock

