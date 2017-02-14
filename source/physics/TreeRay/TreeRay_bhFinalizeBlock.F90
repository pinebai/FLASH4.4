!!****f* source/physics/TreeRay/TreeRay_bhFinalizeBlock
!!
!! NAME
!!
!!  TreeRay_bhFinalizeBlock
!!
!!
!! SYNOPSIS
!!
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!
!! RESULT
!!
!!
!!***

subroutine TreeRay_bhFinalizeBlock(blockno, blkLimits, solnData)
  implicit none
#include "constants.h"
  integer, intent(IN) :: blockno
  integer, dimension(2,MDIM), intent(IN)   :: blkLimits
  real, DIMENSION(:,:,:,:), POINTER :: solnData

end subroutine TreeRay_bhFinalizeBlock
