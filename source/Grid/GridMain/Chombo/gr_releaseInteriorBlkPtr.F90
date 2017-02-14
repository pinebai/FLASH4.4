!!****if* source/Grid/GridMain/Chombo/gr_releaseInteriorBlkPtr
!!
!! NAME
!!  gr_releaseInteriorBlkPtr
!!
!! SYNOPSIS
!!
!!  gr_releaseInteriorBlkPtr(integer(in) :: blockID
!!                           real(pointer) :: dataPtr(:,:,:,:)
!!                           integer(in) :: gridDataStruct)
!!  
!! DESCRIPTION 
!!  Releases a pointer to a block.
!!  
!! ARGUMENTS 
!!
!!  blockID - The block being pointed to.
!!  dataPtr - Pointer to be released.
!!  gridDataStruct - The type of block being pointed to.
!!
!!
!! SEE ALSO
!!  gr_getInteriorBlkPtr
!!***

!!REORDER(4): tmpPtr, dataPtr
#include "constants.h"
#include "Flash.h"

subroutine gr_releaseInteriorBlkPtr(blockID,dataPtr,gridDataStruct)
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
       Grid_getBlkIndexLimits, Grid_getNumVars
  use Driver_interface, ONLY : Driver_abortFlash
  implicit none
  integer,intent(in) :: blockID
  real, pointer :: dataPtr(:,:,:,:)
  integer, intent(in) :: gridDataStruct
  real, dimension(:,:,:,:), pointer :: tmpPtr
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  integer :: nVar

  call Grid_getNumVars(gridDataStruct,nVar)
  call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)
  call Grid_getBlkPtr(blockID, tmpPtr, gridDataStruct)
  tmpPtr( &
       1:nVar, &
       blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
       blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS), &
       blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) = dataPtr
  call Grid_releaseBlkPtr(blockID, tmpPtr, gridDataStruct)

  deallocate(dataPtr)
  nullify(dataPtr)

end subroutine gr_releaseInteriorBlkPtr
