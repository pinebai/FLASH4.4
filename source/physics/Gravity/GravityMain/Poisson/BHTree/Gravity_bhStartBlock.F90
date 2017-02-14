!!****if* source/physics/Gravity/GravityMain/Poisson/BHTree/Gravity_bhStartBlock
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
  use Gravity_data, ONLY : useGravity, grv_bhMeanBlockAccErrInv
  implicit none
#include "constants.h"
#include "Flash.h"
#include "FortranLangFeatures.fh"
  integer, intent(IN) :: blockno
  integer, dimension(2,MDIM), intent(IN)   :: blkLimits
  real, DIMENSION(:,:,:,:), POINTER_INTENT_IN :: solnData
  integer :: ii, jj, kk

  if (.not. useGravity) return

  grv_bhMeanBlockAccErrInv = 0.0
  do kk = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS) ! this block
    do jj = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
      do ii = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
        grv_bhMeanBlockAccErrInv = grv_bhMeanBlockAccErrInv &
        & + solnData(ACEI_VAR,ii,jj,kk)
      enddo
    enddo
  enddo
  grv_bhMeanBlockAccErrInv = grv_bhMeanBlockAccErrInv / ( &
  & (blkLimits(HIGH,IAXIS) - blkLimits(LOW,IAXIS)) * &
  & (blkLimits(HIGH,JAXIS) - blkLimits(LOW,JAXIS)) * &
  & (blkLimits(HIGH,KAXIS) - blkLimits(LOW,KAXIS)) )



  return
end subroutine Gravity_bhStartBlock

