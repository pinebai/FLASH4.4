!!****f* source/physics/TreeRay/TreeRay_bhCheckAccuracy
!!
!! NAME
!!
!!  TreeRay_bhCheckAccuracy
!!
!!
!! SYNOPSIS
!!
!!  logical res = TreeRay_bhCheckAccuracy(
!!                          integer(in)    :: blockno,
!!                          integer(in)    :: point(MDIM),
!!                          integer(in)    :: blkLimits(2,MDIM),
!!                          real,pointer   :: solnData(:,:,:,:)
!!        )
!!
!! DESCRIPTION
!!
!!  Evaluates quantities calculated by TreeRay within a given block 
!!  are accurate enough and do not need to be recalculated.
!!
!! ARGUMENTS
!!
!!  blockno     : number of block into which the target cell belongs
!!  point       : indeces of the target cell in the block
!!  blkLimits   : limits of indeces in the block
!!  solnData    : solution data from the grid
!!
!! RESULT
!!
!!  Returns TRUE if quantities calculated by TreeRay in a given block 
!!  are accurate enough and do not need to be recalculated.
!!
!!***

logical function TreeRay_bhCheckAccuracy(block, point, blkLimits, solnData)
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
  integer, intent(IN) :: block
  integer, dimension(MDIM), intent(IN) :: point
  integer, dimension(2,MDIM), intent(IN)   :: blkLimits
  real, DIMENSION(:,:,:,:), POINTER :: solnData

  TreeRay_bhCheckAccuracy = .true.
  return
end function TreeRay_bhCheckAccuracy


