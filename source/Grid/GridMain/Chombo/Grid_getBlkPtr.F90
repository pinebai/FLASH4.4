!!****if* source/Grid/GridMain/Chombo/Grid_getBlkPtr
!!
!! NAME
!!  Grid_getBlkPtr
!!
!! SYNOPSIS
!!
!!  Grid_getBlkPtr(integer(IN)            :: blockID,
!!                 real(pointer)(:,:,:,:) :: dataPtr,
!!                 integer(IN),optional   :: gridDataStruct)
!!  
!! DESCRIPTION 
!!  
!!  Gets a pointer to a single block of simulation data from the
!!  specified Grid data structure. The block includes guard cells.
!!  If the optional argument "gridDataStructure" is not specified,
!!  it returns a block from cell centered data structure.
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
!!                   SCRATCH space for saving previous timestep data if needed
!!
!!
!!
!! NOTES
!!
!!  Grid_getBlkPtr is an accessor function that passes a pointer
!!  as an argument and requires an explicit interface for most compilers.
!!
!!  Don't forget to call Grid_releaseBlkPtr when you are finished with it!
!!
!!***

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

#include "constants.h"
#include "Flash.h"

subroutine Grid_getBlkPtr(blockID, dataPtr, gridDataStruct)
  use iso_c_binding, only : c_ptr, c_int, c_f_pointer, c_associated
  use flash_ftypes, ONLY : box_info_t
  use chombo_f_c_interface, ONLY : ch_get_blk_ptr, ch_get_box_info
  use Driver_interface, ONLY : Driver_abortFlash
  implicit none
  integer, intent(in) :: blockID

  !dataPtr should really have type real(c_double), but since we know we are
  !promoting FLASH reals to double precision this is OK.
  real, dimension(:,:,:,:), pointer :: dataPtr
  integer, optional,intent(in) :: gridDataStruct

  integer, dimension(MDIM) :: blkSizeInclGC
  integer(c_int), dimension(MDIM) :: loIndex, hiIndex, gCells
  integer(c_int) :: blkID, gds
  integer :: nVar
  type(c_ptr) :: dataMemoryAddress
  logical :: validGridDataStruct
  type(box_info_t) :: boxInfo

#ifdef DEBUG_GRID
  if(present(gridDataStruct)) then
     validGridDataStruct = .false.
     validGridDataStruct= (gridDataStruct == CENTER).or.validGridDataStruct
     validGridDataStruct= (gridDataStruct == FACEX).or.validGridDataStruct
     validGridDataStruct= (gridDataStruct == FACEY).or.validGridDataStruct
     validGridDataStruct= (gridDataStruct == FACEZ).or.validGridDataStruct
     validGridDataStruct= (gridDataStruct == SCRATCH).or.validGridDataStruct
     validGridDataStruct= (gridDataStruct == SCRATCH_CTR).or.validGridDataStruct
     validGridDataStruct= (gridDataStruct == SCRATCH_FACEX).or.validGridDataStruct
     validGridDataStruct= (gridDataStruct == SCRATCH_FACEY).or.validGridDataStruct
     validGridDataStruct= (gridDataStruct == SCRATCH_FACEZ).or.validGridDataStruct
#ifdef FLASH_GRID_PARAMESH
     validGridDataStruct= (gridDataStruct == WORK).or.validGridDataStruct
#endif
     
     if(.not.validGridDataStruct) then
        print *, "Grid_getBlkPtr: gridDataStruct set to improper value"
        print *, "gridDataStruct must = CENTER,FACEX,FACEY,FACEZ," // &
             "WORK or SCRATCH (defined in constants.h)"
        call Driver_abortFlash("gridDataStruct must be one of CENTER,FACEX,FACEY,FACEZ,SCRATCH (see constants.h)")
     end if
  end if
#endif


  !Convert to c_int type.
  blkID = blockID
  if(present(gridDataStruct)) then
     gds = gridDataStruct
  else
     gds = CENTER
  end if


  select case (gds)
  case (CENTER)
     nVar = NUNK_VARS
  case (FACEX, FACEY, FACEZ)
     nVar = NFACE_VARS
  case (SCRATCH)
     nVar = NSCRATCH_GRID_VARS
  case (SCRATCH_CTR)
     nVar = NSCRATCH_CENTER_VARS
  case (SCRATCH_FACEX)
     nVar = NSCRATCH_FACEX_VARS
  case (SCRATCH_FACEY)
     nVar = NSCRATCH_FACEY_VARS
  case (SCRATCH_FACEZ)
     nVar = NSCRATCH_FACEZ_VARS
  case DEFAULT
     call Driver_abortFlash("Grid data structure not yet handled")
  end select


  if (gds == CENTER .or. gds == SCRATCH_CTR) then
     !Obtain a block's raw memory address.
     call ch_get_blk_ptr(blkID, gds, dataMemoryAddress);
     if (c_associated(dataMemoryAddress).eqv..true.) then

        call ch_get_box_info(blkID, gds, boxInfo);
        hiIndex = boxInfo % highLimits
        loIndex = boxInfo % lowLimits
        gcells = boxInfo % guardcells
        blkSizeInclGC = (hiIndex - loIndex) + (2*gCells) + 1


        !Construct a Fortran pointer object using raw memory address and shape argument.
        call c_f_pointer(dataMemoryAddress, dataPtr, &
             (/blkSizeInclGC(1),blkSizeInclGC(2),blkSizeInclGC(3),nVar/))
     else
        nullify(dataPtr)
     end if

     if (.not.associated(dataPtr)) then
        call Driver_abortFlash("[Grid_getBlkPtr]: Cannot construct Fortran pointer")
     end if
  else
     call Driver_abortFlash("[Grid_getBlkPtr]: Grid structure not supported")
  end if

end subroutine Grid_getBlkPtr
