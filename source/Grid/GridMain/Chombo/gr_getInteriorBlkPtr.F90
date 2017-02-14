!!****if* source/Grid/GridMain/Chombo/gr_getInteriorBlkPtr
!!
!! NAME
!!  gr_getInteriorBlkPtr
!!
!! SYNOPSIS
!!
!!  gr_getInteriorBlkPtr(integer(IN)      :: blockID,
!!                 real(pointer)(:,:,:,:) :: dataPtr,
!!                 integer(IN),optional   :: gridDataStruct)
!!  
!! DESCRIPTION 
!!  
!!  Gets a pointer to a single block of simulation data from the
!!  paramesh data structures when there are no permanent guard cells
!!
!!  When using Paramesh 4 in NO_PERMANENT_GUARDCELLS mode, it is important to
!!  release the block pointer for a block before getting it for another block.
!!  For example if pointer to block 1 is not yet released and the user
!!  tries to get a pointer to block 2, the routine will abort.
!!
!! ARGUMENTS 
!!
!!  blockID : the local blockid
!!
!!  dataPtr : Pointer to the data block
!!
!!  gridDataStruct : optional integer value specifying data structure. 
!!                   The options are defined in constants.h and they are :
!!                   CENTER cell centered variables (default)
!!                   FACEX  face centered variable on faces along IAXIS
!!                   FACEY  face centered variable on faces along JAXIS
!!                   FACEZ  face centered variable on faces along IAXIS
!!
!!
!!
!! NOTES
!!
!!  gr_getInteriorBlkPtr is an accessor function that passes a pointer
!!  as an argument and requires an explicit interface for most compilers.
!!
!!  Don't forget to call Grid_releaseBlkPtr when you are finished with it!
!!
!!***

!!REORDER(4): tmpPtr, dataPtr
#include "constants.h"
#include "Flash.h"

subroutine gr_getInteriorBlkPtr(blockID, dataPtr, gridDataStruct)
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
       Grid_getBlkIndexLimits, Grid_getNumVars
  use Driver_interface, ONLY : Driver_abortFlash
  implicit none
  integer, intent(in) :: blockID
  real, dimension(:,:,:,:), pointer :: dataPtr
  integer, intent(in) :: gridDataStruct

  real, dimension(:,:,:,:), pointer :: tmpPtr
  integer, dimension(MDIM) :: sizeMinusGC
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  integer :: nVar

  call Grid_getNumVars(gridDataStruct,nVar)
  call Grid_getBlkIndexLimits(blockId, blkLimits, blkLimitsGC)
  sizeMinusGC(1:MDIM) = blkLimits(HIGH,1:MDIM) - blkLimits(LOW,1:MDIM) + 1
  if (any(SizeMinusGC(1:MDIM) < 1) .or. (nVar < 1)) then
     call Driver_abortFlash("[gr_getInteriorBlkPtr]: Invalid size")
  end if

  !Important that the user calls gr_releaseInteriorBlkPtr.
  allocate(dataPtr( &
       nVar, &
       sizeMinusGC(IAXIS), &
       sizeMinusGC(JAXIS), &
       sizeMinusGC(KAXIS)))

  call Grid_getBlkPtr(blockID, tmpPtr, gridDataStruct)
  dataPtr = tmpPtr( &
       1:nVar, &
       blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
       blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS), &
       blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))
  call Grid_releaseBlkPtr(blockID, tmpPtr, gridDataStruct)

end subroutine gr_getInteriorBlkPtr
