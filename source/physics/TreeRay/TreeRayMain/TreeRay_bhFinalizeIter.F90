!!****if* source/physics/TreeRay/TreeRayMain/TreeRay_bhFinalizeIter
!!
!! NAME
!!
!!  TreeRay_bhFinalizeIter
!!  
!! SYNOPSIS
!!
!!  TreeRay_bhFinalizeIter()
!!
!! DESCRIPTION
!!
!!
!!***
subroutine TreeRay_bhFinalizeIter()
  use TreeRay_data, ONLY : tr_BhUseTreeRay
  use Grid_interface, ONLY : Grid_getListOfBlocks, &
    Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_releaseBlkPtr, &
    Grid_getDeltas
  use tr_bhLocalInterface, ONLY : tr_bhRadToGas
  use Driver_interface, ONLY : Driver_getDt
  use Eos_interface, ONLY : Eos_wrapped


  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  integer :: i, j, k, bl, blockCount, blockno
  integer :: blockList(MAXBLOCKS)  
  integer, dimension(2,MDIM)   :: blkLimits, blkLimitsGC  
  real :: vol_poc, area_poc, del(MDIM), dt
  real, POINTER, DIMENSION(:,:,:,:) :: solnData
  real, dimension(:), pointer :: solnPoint

  if (.not. tr_bhUseTreeRay) return

  call Driver_getDt(dt)
  call Grid_getListOfBlocks(LEAF,blockList,blockCount)

  do bl = 1, blockCount
    blockno = blockList(bl)

    ! get information about the block
    call Grid_getBlkIndexLimits(blockno,blkLimits,blkLimitsGC)
    call Grid_getDeltas(blockno,del)
    call Grid_getBlkPtr(blockno,solnData,CENTER) ! pointer to the density and gp field


    vol_poc = del(IAXIS)*del(JAXIS)*del(KAXIS)
    area_poc = del(IAXIS)*del(IAXIS)

    do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
      do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
        do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)

          ! sum up contributions from all rays
          !NUV_poc(:) = NUV_poc(:) + eflux(:,ipix) * area_poc
        
          solnPoint => solnData(:,i,j,k)

          call tr_bhRadToGas(solnPoint, vol_poc, area_poc, dt)
        enddo
      enddo
    enddo
    call Eos_wrapped(MODE_DENS_EI, blkLimits, blockno)
    call Grid_releaseBlkPtr(blockno,solnData,CENTER) ! pointer to the density and gp field

  enddo



  return

end subroutine TreeRay_bhFinalizeIter
