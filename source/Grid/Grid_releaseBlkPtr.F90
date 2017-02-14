!!****f* source/Grid/Grid_releaseBlkPtr
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
!!  gridDataStruct - see Grid_getBlkPtr
!!
!! NOTES
!!
!!  Some implementations actually do more than just releasing the pointer.
!!
!! SEE ALSO
!!  Grid_getBlkPtr
!!***

subroutine Grid_releaseBlkPtr(blockId, dataPtr, gridDataStruct)

  implicit none
  integer,intent(in) :: blockId
  real, pointer :: dataPtr(:,:,:,:)
  integer,optional, intent(in) :: gridDataStruct

end subroutine Grid_releaseBlkPtr
