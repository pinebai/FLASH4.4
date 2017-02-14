!!****if* source/Grid/GridMain/Chombo/Grid_releaseBlkPtr
!!
!! NAME
!!  Grid_releaseBlkPtr
!!
!! SYNOPSIS
!!
!!  Grid_releaseBlkPtr(integer(IN)   :: blockId,
!!                     real(pointer) :: dataPtr(:,:,:,:),
!!                     integer(IN),optional :: gridDataStruct)
!!  
!! DESCRIPTION 
!!  Releases a pointer to a block.
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
!!  Grid_getBlkPtr
!!***

subroutine Grid_releaseBlkPtr(blockId, dataPtr, gridDataStruct)

  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

  integer,intent(in) :: blockId
  real, pointer :: dataPtr(:,:,:,:)
  integer,optional, intent(in) :: gridDataStruct

  nullify(dataPtr)

end subroutine Grid_releaseBlkPtr
