!!****if* source/Grid/GridMain/paramesh/paramesh4/gr_getInteriorBlkPtr
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

#include "constants.h"
#include "Flash.h"

  use physicaldata, ONLY : unk, facevarx, facevary, facevarz
  use Driver_interface, ONLY : Driver_abortFlash

#ifdef FLASH_GRID_PARAMESH
  use workspace, ONLY : work
#endif 

  implicit none
  integer, intent(in) :: blockID
  real, dimension(:,:,:,:), pointer :: dataPtr
  integer, intent(in) :: gridDataStruct

  logical :: validGridDataStruct

  integer :: idest, iopt, nlayers, icoord
  logical :: lcc, lfc, lec, lnc, l_srl_only, ldiag
  logical,dimension(NUNK_VARS) :: save_ccMask

#ifndef FL_NON_PERMANENT_GUARDCELLS
  call Driver_abortFlash("interior of blocks can be got only in non permanent gc mode")
#endif

#ifndef FLASH_GRID_PARAMESH
  call Driver_abortFlash("interior of blocks can be got only from Paramesh")
#endif



#ifdef DEBUG_GRID
  validGridDataStruct = .false.
  validGridDataStruct= (gridDataStruct == CENTER).or.validGridDataStruct
  validGridDataStruct= (gridDataStruct == FACEX).or.validGridDataStruct
  validGridDataStruct= (gridDataStruct == FACEY).or.validGridDataStruct
  validGridDataStruct= (gridDataStruct == FACEZ).or.validGridDataStruct
  validGridDataStruct= (gridDataStruct == WORK).or.validGridDataStruct
  
  if(.not.validGridDataStruct) then
     print *, "gr_getInteriorBlkPtr: gridDataStruct set to improper value"
     print *, "gridDataStruct must = CENTER,FACEX,FACEY,FACEZ, or " // &
          "WORK (defined in constants.h)"
     call Driver_abortFlash("gridDataStruct must be one of CENTER,FACEX,FACEY,FACEZ (see constants.h)")
  end if
#endif

  if((blockid<1).or.(blockid>MAXBLOCKS)) then
     print *, 'gr_getInteriorBlkPtr:  invalid blockid ',blockid
     call Driver_abortFlash("[gr_getInteriorBlkPtr] invalid blockid ")
  end if
  
  select case(gridDataStruct)
  case(CENTER)
     dataPtr => unk(:,:,:,:,blockid)
  case(FACEX)
     dataPtr => facevarx(:,:,:,:,blockid)
  case(FACEY)
     dataPtr => facevary(:,:,:,:,blockid)
  case(FACEZ)
     dataPtr => facevarz(:,:,:,:,blockid)
#ifdef FLASH_GRID_PARAMESH
  case(WORK)
     call Driver_abortFlash( &
          "[gr_getInteriorBlkPtr] work array cannot be got as pointer - NOT IMPLEMENTED")
#endif
  end select
  return
end subroutine gr_getInteriorBlkPtr








