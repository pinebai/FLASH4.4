!!****if* source/physics/TreeRay/TreeRayMain/TreeRay_bhStartBlock
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
  use TreeRay_data, ONLY : tr_bhMassRays, tr_bhVolRays, tr_bhSrcfRays, &
    tr_bhEradRays, tr_bhUseTreeRay
  !use tr_osInterface, ONLY : tr_osStartBlock
  use tr_odInterface, ONLY : tr_odStartBlock
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
  integer, intent(IN) :: blockno
  integer, dimension(2,MDIM)   :: blkLimits
  real, DIMENSION(:,:,:,:), POINTER, intent(OUT) :: solnData

  if (.not. tr_bhUseTreeRay) return

  call tr_odStartBlock(blockno, blkLimits, solnData)

  return
end subroutine TreeRay_bhStartBlock

