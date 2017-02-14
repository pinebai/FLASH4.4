!!****if* source/physics/TreeRay/TreeRayMain/TreeRay_bhAccBotNode
!!
!! NAME
!!
!!  TreeRay_bhAccBotNode
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

subroutine TreeRay_bhAccBotNode(blockno, point, blkLimits, locCoords, solnData, &
   & botnode, accnode)
  !use tr_osInterface, ONLY : tr_osAccBotNode
  use tr_odInterface, ONLY : tr_odAccBotNode
  use TreeRay_data, ONLY : tr_bhUseTreeRay
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

  if (.not. tr_bhUseTreeRay) return

  !call tr_osAccBotNode(blockno, point, blkLimits, locCoords, solnData, &
  ! & botnode, accnode)
  call tr_odAccBotNode(blockno, point, blkLimits, locCoords, solnData, &
   & botnode, accnode)

end subroutine TreeRay_bhAccBotNode
