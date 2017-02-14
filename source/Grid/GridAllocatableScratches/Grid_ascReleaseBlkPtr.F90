!!****if* source/Grid/GridAllocatableScratches/Grid_ascReleaseBlkPtr
!!
!! NAME
!!  Grid_ascReleaseBlkPtr
!!
!! SYNOPSIS
!!
!!  call Grid_ascReleaseBlkPtr(integer(IN)   :: blockId,
!!                     real(pointer) :: dataPtr(:,:,:,:),
!!                     integer(IN),optional :: gridDataStruct)
!!  
!! DESCRIPTION 
!!  Releases a pointer to allocatable scratch data for a block.
!!  
!! ARGUMENTS 
!!
!!  blockId - ID of the block, should be the same ID was used in the
!!            corresponding Grid_getBlkPtr call.
!!  dataPtr - Pointer to be released.
!!  gridDataStruct - an optional argument that designates the type of grid data
!!                   structure to handle (i.e. facevar, unknown, scratch...)
!!
!!
!! NOTES
!!
!!  This implementation actually does more than just releasing the pointer.
!!
!! SEE ALSO
!!  Grid_ascGetBlkPtr
!!  Grid_releaseBlkPtr
!!***

subroutine Grid_ascReleaseBlkPtr(blockId, dataPtr, gridDataStruct)

#include "Flash.h"
#ifdef Grid_ascReleaseBlkPtr
! disabling drift macro expansion, see: drift
#undef Grid_ascReleaseBlkPtr
#endif


#include "constants.h"


  implicit none

  integer,intent(in) :: blockId
  real, pointer :: dataPtr(:,:,:,:)
  integer,optional, intent(in) :: gridDataStruct

  integer :: gds

  if(present(gridDataStruct)) then
     gds = gridDataStruct
  else
     gds = CENTER
  end if


  ! always destroy the pointer, because the other users will have their 
  ! own
  nullify(dataPtr)

end subroutine Grid_ascReleaseBlkPtr

subroutine Grid_ascReleaseBlk5Ptr(blockId, data5Ptr, gridDataStruct)

#include "FortranLangFeatures.fh"

  implicit none

  integer,intent(in) :: blockId
  real, POINTER_INTENT_OUT :: data5Ptr(:,:,:,:,:)
  integer,optional, intent(in) :: gridDataStruct

  integer :: gds

  if(present(gridDataStruct)) then
     gds = gridDataStruct
  else
     gds = CENTER
  end if

  nullify(data5Ptr)

end subroutine Grid_ascReleaseBlk5Ptr
