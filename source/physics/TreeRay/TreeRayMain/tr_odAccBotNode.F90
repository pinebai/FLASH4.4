!!****if* source/physics/TreeRay/TreeRayMain/tr_odAccBotNode
!!
!! NAME
!!
!!  tr_odAccBotNode
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

subroutine tr_odAccBotNode(blockno, point, blkLimits, locCoords, solnData, &
   & botnode, accnode)
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  integer, intent(IN) :: blockno
  integer, dimension(MDIM), intent(IN) :: point
  integer, dimension(2,MDIM)   :: blkLimits
  real, dimension(:,:,:), intent(IN) :: locCoords
  real, DIMENSION(:,:,:,:), POINTER, intent(IN) :: solnData
  real, dimension(:), intent(IN) :: botnode
  real, dimension(:), intent(INOUT) :: accnode

end subroutine tr_odAccBotNode
