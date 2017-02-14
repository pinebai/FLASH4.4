!!****if* source/physics/Gravity/GravityMain/Poisson/BHTree/Gravity_bhFinalizeBlock
!!
!! NAME
!!
!!  Gravity_bhFinalizeBlock
!!
!!
!! SYNOPSIS
!!
!!   call Gravity_bhFinalizeBlock(
!!                          integer(in)    :: blockno,
!!                          integer(in)    :: blkLimits(2,MDIM),
!!                          real,pointer   :: solnData(:,:,:,:)
!!        )
!!
!! DESCRIPTION
!!
!!  Called during the Tree Walk after the tree is traversed for all
!!  cells of the block. In this version empty - no need to do 
!!  anything here for the calculation of the gravitational potential.
!!
!! ARGUMENTS
!!
!!  blockno     : number of block for which the tree walk was just finished
!!  blkLimits   : limits of indeces in the block
!!  solnData    : solution data from the grid
!!
!!
!!***

subroutine Gravity_bhFinalizeBlock(blockno, blkLimits, solnData)
  use Gravity_data, ONLY : useGravity
  implicit none
#include "constants.h"
  integer, intent(IN) :: blockno
  integer, dimension(2,MDIM), intent(IN)   :: blkLimits
  real, DIMENSION(:,:,:,:), POINTER :: solnData

  if (.not. useGravity) return

end subroutine Gravity_bhFinalizeBlock
