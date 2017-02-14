!!****if* source/Grid/localAPI/gr_getInteriorBlkPtr
!!
!! NAME
!!  gr_getInteriorBlkPtr
!!
!! SYNOPSIS
!!
!!  gr_getInteriorBlkPtr(integer(IN)            :: blockID,
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

!!REORDER(5): unk, facevar[xyz]
!!FOR FUTURE: Add REORDER for unk, facevar[xyz]1, etc.?

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

subroutine gr_getInteriorBlkPtr(blockID,dataPtr, gridDataStruct)

  implicit none
  integer, intent(in) :: blockID
  real, dimension(:,:,:,:), pointer :: dataPtr
  integer, intent(in) :: gridDataStruct
  
  return
end subroutine gr_getInteriorBlkPtr








