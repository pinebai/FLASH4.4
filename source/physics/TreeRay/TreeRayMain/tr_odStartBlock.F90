!!****if* source/physics/TreeRay/TreeRayMain/tr_odStartBlock
!!
!! NAME
!!
!!  tr_odStartBlock
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

subroutine tr_odStartBlock(blockno, blkLimits, solnData)
  implicit none
#include "constants.h"
#include "Flash.h"
  integer, intent(IN) :: blockno
  integer, dimension(2,MDIM)   :: blkLimits
  real, DIMENSION(:,:,:,:), POINTER, intent(OUT) :: solnData

  return
end subroutine tr_odStartBlock

