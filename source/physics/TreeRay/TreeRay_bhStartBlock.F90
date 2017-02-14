!!****f* source/physics/TreeRay/TreeRay_bhStartBlock
!!
!! NAME
!!
!!  TreeRay_bhStartBlock
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

subroutine TreeRay_bhStartBlock(blockno, blkLimits, solnData)
  implicit none
#include "constants.h"
#include "Flash.h"
  integer, intent(IN) :: blockno
  integer, dimension(2,MDIM), intent(IN)   :: blkLimits
  real, DIMENSION(:,:,:,:), POINTER :: solnData

  return
end subroutine TreeRay_bhStartBlock

