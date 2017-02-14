!!****if* source/Grid/GridMain/paramesh/paramesh4/gr_releaseInteriorBlkPtr
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

subroutine gr_releaseInteriorBlkPtr(blockID,dataPtr,gridDataStruct)

  implicit none
  integer,intent(in) :: blockID
  real, pointer :: dataPtr(:,:,:,:)
  integer, intent(in) :: gridDataStruct

  nullify(dataPtr)

end subroutine gr_releaseInteriorBlkPtr
