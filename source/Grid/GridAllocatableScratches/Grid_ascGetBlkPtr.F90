!!****if* source/Grid/GridAllocatableScratches/Grid_ascGetBlkPtr
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

!!REORDER(5):gr_ascArrayScratch, gr_ascArrayScratch_ctr, gr_ascArrayScratch_facevar[xyz]
!!REORDER(4): dataPtr

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

subroutine Grid_ascGetBlkPtr(blockID,dataPtr, gridDataStruct)

#include "constants.h"
#include "Flash.h"

  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_ascModule, ONLY : gr_ascArrayScratch,gr_ascArrayCenter,gr_ascArrayScratch_ctr,&
       gr_ascArrayScratch_facevarx,gr_ascArrayScratch_facevary,gr_ascArrayScratch_facevarz, &
       gr_ascArrayFacevarx,gr_ascArrayFacevary,gr_ascArrayFacevarz, &
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
     validGridDataStruct= (gridDataStruct == FACEX).or.validGridDataStruct
     validGridDataStruct= (gridDataStruct == FACEY).or.validGridDataStruct
     validGridDataStruct= (gridDataStruct == FACEZ).or.validGridDataStruct
#ifdef FLASH_GRID_PARAMESH
     validGridDataStruct= (gridDataStruct == WORK).or.validGridDataStruct
#endif
     
     if(.not.validGridDataStruct) then
        print *, "Grid_ascGetBlkPtr: gridDataStruct set to improper value"
        print *, "gridDataStruct must = SCRATCH_CTR,SCRATCH_FACEX,_FACEY,_FACEZ," // &
             "or SCRATCH, etc. (defined in constants.h)"
        call Driver_abortFlash("gridDataStruct must be SCRATCH or SCRATCH_*, etc. (see constants.h)")
     end if
  end if
  if((blockid<1).or.(blockid>MAXBLOCKS)) then
     print *, 'Grid_ascGetBlkPtr:  invalid blockid ',blockid
     call Driver_abortFlash("[Grid_ascGetBlkPtr] invalid blockid ")
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
     print*, "Grid_ascGetBlkPtr: blockID",blockID," -> invalid leafNo",leafNo
     call Driver_abortFlash("Grid_ascGetBlkPtr: invalid leafNo!")
  end if
#endif

#include "FortranLangFeatures.fh"

#define SETDATAPTR(arr) medPtr=>arr; call AssoMed(dataPtr,medPtr,leafNo)

  select case (gds)
  case(SCRATCH)
     SETDATAPTR(gr_ascArrayScratch)
  case(CENTER)
     SETDATAPTR(gr_ascArraycenter)           
  case(SCRATCH_CTR)
     SETDATAPTR(gr_ascArrayScratch_ctr)           
  case(SCRATCH_FACEX)
     SETDATAPTR(gr_ascArrayScratch_facevarx)
  case(SCRATCH_FACEY)
     SETDATAPTR(gr_ascArrayScratch_facevary)           
  case(SCRATCH_FACEZ)
     SETDATAPTR(gr_ascArrayScratch_facevarz)           
  case(FACEX)
     SETDATAPTR(gr_ascArrayFacevarx)
  case(FACEY)
     SETDATAPTR(gr_ascArrayFacevary)           
  case(FACEZ)
     SETDATAPTR(gr_ascArrayFacevarz)           
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

end subroutine Grid_ascGetBlkPtr


subroutine Grid_ascGetBlk5Ptr(blockID,data5Ptr, gridDataStruct)

#include "constants.h"
#include "Flash.h"

  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_ascModule, ONLY : gr_ascArray5Center, &
       gr_ascArray5Scratch_ctr, &
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
        print *, "Grid_ascGetBlk5Ptr: gridDataStruct set to improper value"
        print *, "gridDataStruct must be CENTER or SCRATCH_CTR," // &
             "or SCRATCH (defined in constants.h)"
        call Driver_abortFlash("gridDataStruct must be CENTER or SCRATCH_CTR(see constants.h)")
     end if
  end if
  if((blockid<1).or.(blockid>MAXBLOCKS)) then
     print *, 'Grid_ascGetBlk5Ptr:  invalid blockid ',blockid
     call Driver_abortFlash("[Grid_ascGetBlk5Ptr] invalid blockid ")
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
     print*, "Grid_ascGetBlk5Ptr: blockID",blockID," -> invalid leafNo",leafNo
     call Driver_abortFlash("Grid_ascGetBlk5Ptr: invalid leafNo!")
  end if
#endif


#define SETDATA5PTR(arr) med5Ptr=>arr; call AssoMed(data5Ptr,med5Ptr,leafNo)

  select case (gds)
  case(CENTER)
     SETDATA5PTR(gr_ascArray5Center)           
  case(SCRATCH_CTR)
     SETDATA5PTR(gr_ascArray5Scratch_ctr)           
  case DEFAULT
     print *, 'TRIED TO GET SOMETHING OTHER THAN CENTER OR SCRATCH_CTR. NOT VALID.'
     call Driver_abortFlash('Grid_ascGetBlk5Ptr: TRIED TO GET SOMETHING OTHER THAN CENTER OR SCRATCH_CTR. NOT VALID.')
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

end subroutine Grid_ascGetBlk5Ptr
