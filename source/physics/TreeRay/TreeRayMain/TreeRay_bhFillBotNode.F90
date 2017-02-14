!!****if* source/physics/TreeRay/TreeRayMain/TreeRay_bhFillBotNode
!!
!! NAME
!!
!!  TreeRay_bhFillBotNode
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

subroutine TreeRay_bhFillBotNode(blockno, point, blkLimits, solnData, botnode)
  use TreeRay_data, ONLY : tr_BhUseTreeRay
  !use tr_osInterface, ONLY : tr_osFillBotNode
  use tr_odInterface, ONLY : tr_odFillBotNode
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  integer, intent(IN) :: blockno
  integer, dimension(MDIM), intent(IN) :: point
  integer, dimension(2,MDIM)   :: blkLimits
  real, DIMENSION(:,:,:,:), POINTER, intent(IN) :: solnData
  real, dimension(:), intent(OUT) :: botnode
  real :: dvol


  if (.not. tr_bhUseTreeRay) return

  !call tr_osFillBotNode(blockno, point, blkLimits, solnData, botnode)
  call tr_odFillBotNode(blockno, point, blkLimits, solnData, botnode)

  return
end subroutine TreeRay_bhFillBotNode

