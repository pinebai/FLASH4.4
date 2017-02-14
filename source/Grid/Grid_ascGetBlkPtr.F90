!!****f* source/Grid/Grid_ascGetBlkPtr
!!
!! NAME
!!  Grid_ascGetBlkPtr
!!
!! SYNOPSIS
!!
!!  call Grid_ascGetBlkPtr(integer(IN)            :: blockID,
!!                 real(pointer)(:,:,:,:) :: dataPtr,
!!                 integer(IN),optional   :: gridDataStruct)
!!  
!! DESCRIPTION 
!!  
!!  Gets a pointer to allocatable scratch data for a single block of simulation data from the
!!  specified Grid data structure. The scratch data may include (some or all) guard cells,
!!  this depends on the preceding Grid_ascAllocMem call.
!!  If the optional argument "gridDataStructure" is not specified,
!!  it returns a block from cell centered data structure.
!!
!! ARGUMENTS 
!!
!!  blockID : the local blockid
!!
!!  dataPtr : Pointer to the data block
!!
!!  gridDataStruct : optional integer value specifying data structure. 
!!                   The options are defined in constants.h and they are :
!!                   SCRATCH scratch space that can fit cell and face centered variables
!!                   SCRATCH_CTR scratch space for cell centered variables
!!                   SCRATCH_FACEX scratch space facex variables
!!                   SCRATCH_FACEY scratch space facey variables
!!                   SCRATCH_FACEZ scratch space facez variables
!!
!!
!!
!! NOTES
!!
!!  Grid_ascGetBlkPtr is an accessor function that passes a pointer
!!  as an argument and requires an explicit interface for most compilers.
!!
!!  Don't forget to call Grid_ascReleaseBlkPtr when you are finished with it!
!!
!! SEE ALSO
!!  Grid_ascAllocMem
!!  Grid_ascReleaseBlkPtr
!!  Grid_getBlkPtr
!!***

subroutine Grid_ascGetBlkPtr(blockID,dataPtr, gridDataStruct)

  implicit none
  integer, intent(in) :: blockID
  real, dimension(:,:,:,:), pointer :: dataPtr
  integer, optional,intent(in) :: gridDataStruct

end subroutine Grid_ascGetBlkPtr


subroutine Grid_ascGetBlk5Ptr(blockID,data5Ptr, gridDataStruct)

  implicit none
  integer, intent(in) :: blockID
  real, dimension(:,:,:,:,:), pointer :: data5Ptr
  integer, optional,intent(in) :: gridDataStruct
end subroutine Grid_ascGetBlk5Ptr
