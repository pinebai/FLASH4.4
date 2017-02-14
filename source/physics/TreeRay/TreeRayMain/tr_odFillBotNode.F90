!!****if* source/physics/TreeRay/TreeRayMain/tr_odFillBotNode
!!
!! NAME
!!
!!  tr_odFillBotNode
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

subroutine tr_odFillBotNode(blockno, point, blkLimits, solnData, botnode)
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
  integer, intent(IN) :: blockno
  integer, dimension(MDIM), intent(IN) :: point
  integer, dimension(2,MDIM)   :: blkLimits
  real, DIMENSION(:,:,:,:), POINTER, intent(IN) :: solnData
  real, dimension(:), intent(OUT) :: botnode

  return
end subroutine tr_odFillBotNode

