!!****if* source/physics/Hydro/HydroMain/hy_memGetBlkPtr
!!
!! NAME
!!  hy_memGetBlkPtr
!!
!! SYNOPSIS
!!
!!  call hy_memGetBlkPtr(integer(IN)            :: blockID,
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
!!  hy_memGetBlkPtr is an accessor function that passes a pointer
!!  as an argument and requires an explicit interface for most compilers.
!!
!!  Don't forget to call hy_memReleaseBlkPtr when you are finished with it!
!!
!!***

!!REORDER(5):hy_memArrayScratch, hy_memArrayScratch_ctr, hy_memArrayScratch_facevar[xyz]
!!REORDER(4): dataPtr

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

subroutine hy_memGetBlkPtr(blockID,dataPtr, gridDataStruct)

#include "constants.h"
#include "Flash.h"

  use Driver_interface, ONLY : Driver_abortFlash
  use hy_memScratchData, ONLY : hy_memArrayScratch,hy_memArrayCenter,hy_memArrayScratch_ctr,&
       hy_memArrayScratch_facevarx,hy_memArrayScratch_facevary,hy_memArrayScratch_facevarz, &
       bIdToLeafNo



  implicit none
  integer, intent(in) :: blockID
  real, dimension(:,:,:,:), pointer :: dataPtr
  integer, optional,intent(in) :: gridDataStruct

  real, dimension(:,:,:,:,:), pointer :: medPtr
  integer :: gds
  integer :: leafNo
#ifdef DEBUG_GRID
  logical :: validGridDataStruct
#endif

#ifdef DEBUG_GRID
  if(present(gridDataStruct)) then
     validGridDataStruct = .false.
     validGridDataStruct= (gridDataStruct == SCRATCH).or.validGridDataStruct
     validGridDataStruct= (gridDataStruct == CENTER).or.validGridDataStruct
     validGridDataStruct= (gridDataStruct == SCRATCH_CTR).or.validGridDataStruct
     validGridDataStruct= (gridDataStruct == SCRATCH_FACEX).or.validGridDataStruct
     validGridDataStruct= (gridDataStruct == SCRATCH_FACEY).or.validGridDataStruct
     validGridDataStruct= (gridDataStruct == SCRATCH_FACEZ).or.validGridDataStruct
#ifdef FLASH_GRID_PARAMESH
     validGridDataStruct= (gridDataStruct == WORK).or.validGridDataStruct
#endif
     
     if(.not.validGridDataStruct) then
        print *, "hy_memGetBlkPtr: gridDataStruct set to improper value"
        print *, "gridDataStruct must = SCRATCH_CTR,SCRATCH_FACEX,_FACEY,_FACEZ," // &
             "or SCRATCH (defined in constants.h)"
        call Driver_abortFlash("gridDataStruct must be SCRATCH or SCRATCH_* (see constants.h)")
     end if
  end if
  if((blockid<1).or.(blockid>MAXBLOCKS)) then
     print *, 'hy_memGetBlkPtr:  invalid blockid ',blockid
     call Driver_abortFlash("[hy_memGetBlkPtr] invalid blockid ")
  end if
#endif

  if(present(gridDataStruct)) then
     gds = gridDataStruct
  else
     gds = CENTER
  end if

  leafNo = bIdToLeafNo(blockID)
#ifdef DEBUG_GRID
  if (.NOT.(leafNo > 0 .AND. leafNo .LE. MAXBLOCKS)) then
     print*, "hy_memGetBlkPtr: blockID",blockID," -> invalid leafNo",leafNo
     call Driver_abortFlash("hy_memGetBlkPtr: invalid leafNo!")
  end if
#endif

#include "FortranLangFeatures.fh"

#define SETDATAPTR(arr) medPtr=>arr; call AssoMed(dataPtr,medPtr,leafNo)

  select case (gds)
  case(SCRATCH)
     SETDATAPTR(hy_memArrayScratch)
  case(CENTER)
     SETDATAPTR(hy_memArraycenter)           
  case(SCRATCH_CTR)
     SETDATAPTR(hy_memArrayScratch_ctr)           
  case(SCRATCH_FACEX)
     SETDATAPTR(hy_memArrayScratch_facevarx)
  case(SCRATCH_FACEY)
     SETDATAPTR(hy_memArrayScratch_facevary)           
  case(SCRATCH_FACEZ)
     SETDATAPTR(hy_memArrayScratch_facevarz)           
  case DEFAULT
     print *, 'TRIED TO GET SOMETHING OTHER THAN SCRATCH OR SCRATCH_CTR OR SCRATCH_FACE[XYZ]. NOT VALID.'
  end select

  return

contains
  subroutine AssoMed(pp, mm, leafNo)
    real,POINTER_INTENT_OUT :: pp(:,:,:,:)
    real,POINTER_INTENT_IN  :: mm(:,:,:,:,:)
    integer,intent(in) :: leafNo
    call AssoFin(pp,mm(:,:,:,:,leafNo),lbound(mm,1),lbound(mm,2),lbound(mm,3),lbound(mm,4))
  end subroutine AssoMed

  subroutine AssoFin(pp, dd, lb1,lb2,lb3,lb4)
    real,POINTER_INTENT_OUT :: pp(:,:,:,:)
    integer, intent(in) :: lb1,lb2,lb3,lb4
    real,   intent(in),target :: dd(lb1:,lb2:,lb3:,lb4:)
    pp => dd
  end subroutine AssoFin

end subroutine hy_memGetBlkPtr


subroutine hy_memGetBlk5Ptr(blockID,data5Ptr, gridDataStruct)

#include "constants.h"
#include "Flash.h"

  use Driver_interface, ONLY : Driver_abortFlash
  use hy_memScratchData, ONLY : hy_memArray5Center, &
       hy_memArray5Scratch_ctr, &
       bIdToLeafNo



  implicit none
  integer, intent(in) :: blockID
  real, dimension(:,:,:,:,:), pointer :: data5Ptr
  integer, optional,intent(in) :: gridDataStruct

  real, dimension(:,:,:,:,:,:), pointer :: med5Ptr
  integer :: gds
  integer :: leafNo
#ifdef DEBUG_GRID
  logical :: validGridDataStruct
#endif

#ifdef DEBUG_GRID
  if(present(gridDataStruct)) then
     validGridDataStruct = .false.
     validGridDataStruct= (gridDataStruct == CENTER).or.validGridDataStruct
     validGridDataStruct= (gridDataStruct == SCRATCH_CTR).or.validGridDataStruct
     
     if(.not.validGridDataStruct) then
        print *, "hy_memGetBlk5Ptr: gridDataStruct set to improper value"
        print *, "gridDataStruct must be CENTER or SCRATCH_CTR," // &
             "or SCRATCH (defined in constants.h)"
        call Driver_abortFlash("gridDataStruct must be CENTER or SCRATCH_CTR(see constants.h)")
     end if
  end if
  if((blockid<1).or.(blockid>MAXBLOCKS)) then
     print *, 'hy_memGetBlk5Ptr:  invalid blockid ',blockid
     call Driver_abortFlash("[hy_memGetBlk5Ptr] invalid blockid ")
  end if
#endif

  if(present(gridDataStruct)) then
     gds = gridDataStruct
  else
     gds = CENTER
  end if

  leafNo = bIdToLeafNo(blockID)
#ifdef DEBUG_GRID
  if (.NOT.(leafNo > 0 .AND. leafNo .LE. MAXBLOCKS)) then
     print*, "hy_memGetBlk5Ptr: blockID",blockID," -> invalid leafNo",leafNo
     call Driver_abortFlash("hy_memGetBlk5Ptr: invalid leafNo!")
  end if
#endif


#define SETDATA5PTR(arr) med5Ptr=>arr; call AssoMed(data5Ptr,med5Ptr,leafNo)

  select case (gds)
  case(CENTER)
     SETDATA5PTR(hy_memArray5Center)           
  case(SCRATCH_CTR)
     SETDATA5PTR(hy_memArray5Scratch_ctr)           
  case DEFAULT
     print *, 'TRIED TO GET SOMETHING OTHER THAN CENTER OR SCRATCH_CTR. NOT VALID.'
     call Driver_abortFlash('hy_memGetBlk5Ptr: TRIED TO GET SOMETHING OTHER THAN CENTER OR SCRATCH_CTR. NOT VALID.')
  end select

  return

contains
  subroutine AssoMed(pp, mm, leafNo)
    real,POINTER_INTENT_OUT :: pp(:,:,:,:,:)
    real,POINTER_INTENT_IN  :: mm(:,:,:,:,:,:)
    integer,intent(in) :: leafNo
    call AssoFin(pp,mm(:,:,:,:,:,leafNo),lbound(mm,1),lbound(mm,2),lbound(mm,3),lbound(mm,4))
  end subroutine AssoMed

  subroutine AssoFin(pp, dd, lb1,lb2,lb3,lb4)
    real,POINTER_INTENT_OUT :: pp(:,:,:,:,:)
    integer, intent(in) :: lb1,lb2,lb3,lb4
    real,   intent(in),target :: dd(lb1:,lb2:,lb3:,lb4:,:)
    pp => dd
  end subroutine AssoFin

end subroutine hy_memGetBlk5Ptr
