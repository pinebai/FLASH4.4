!!****if* source/Grid/GridMain/Grid_GCPutScratch
!!
!! NAME
!!  Grid_GCPutScratch
!!
!! SYNOPSIS
!!
!!  Grid_GCPutScratch(integer(IN) :: gridDataStruct,
!!                    logical(IN) :: needSetup,
!!                    logical(IN) :: releaseSetup,
!!           optional,integer(IN) :: blkCount,
!!           optional,integer(IN) :: blkList(:),
!!           optional,integer(IN) :: indCount,
!!           optional,integer(IN) :: indList(:),
!!           optional,integer(IN) :: gcCnt(:),
!!           optional,logical(IN) :: putData)
!!  
!!  
!! DESCRIPTION
!!  
!!  This routine stores guardcells of the specified indices of 
!!  grid data structures from a list of blocks. If putData is .false.
!!  then the routine is said to be called in initialize only mode.
!!  In this mode "needSetUp" must be true. If putData is missing, it is 
!!  assumed to be true. 
!!  If the argument needSetup is true, then all the optional 
!!  arguments must be included in the call, and the routine will
!!  make a call to gr_allocGCScratch function. The first call to 
!!  Grid_GCPutScratch in a simulation must be made with needSetup=.true.
!!  for correct operation. 
!!  The subsequent calls to Grid_GCPutScratch can be made without the optional
!!  arguments by setting needSetup to false, as long the values of the
!!  optional arguments provided in the call with needSetup=true
!!  are unchanged in the simulation.
!!
!!  If there is any change in the grid, list of blocks, or the number
!!  or order of indices, then both needSetup and releaseSetup must be set
!!  to true, and appropriate values must be provided in the optional 
!!  argument, which must all be present.
!!
!! ARGUMENTS
!!            
!!  gridDataStruct : integer value specifying data structure. 
!!                   The options are defined in constants.h, the ones
!!                   relevant to this routine are :
!!                   CENTER cell centered variables (default)
!!                   FACEX  face centered variable on faces along IAXIS
!!                   FACEY  face centered variable on faces along JAXIS
!!                   FACEZ  face centered variable on faces along IAXIS
!!  needSetup      : if true a call to gr_allocateGuardScratchSpace will be
!!                   made
!!  releaseSetup   : if true a call to gr_releaseGuardScratchSpace will be
!!                   made before the call to gr_allocateGuardScratchSpace
!!  blkCount       : count of blocks whose guard cells are to be saved
!!  blkList        : list of blocks whose guard cells are to be saved
!!  indCount       : count of indices in the data structure to be saved
!!  indList        : list of indices in the data structure to be saved
!!  gcCnt          : the count of guardcells along each dimension to be saved             
!!  putData        : indicates whether the routine is called in intialize
!!                   mode
!!
!!
!!  NOTES
!!   variables that start with "gr_" are variables of Grid unit scope
!!   and are stored in the fortran module Grid_data. Variables are not
!!   starting with gr_ are local variables or arguments passed to the 
!!   routine.
!!
!!  SEE ALSO
!!   
!!   gr_GCAllocScratch and gr_GCTransferOneBlk
!!
!!***

!!REORDER(4): solnData

subroutine Grid_GCPutScratch(gridDataStruct,needSetup,releaseSetup,&
     &blkCnt,blkList,indCnt,indList, gcCnt, putData)

  use gr_GCScratchData
  use Grid_interface, ONLY :  Grid_getBlkIndexLimits, Grid_getBlkPtr,&
       Grid_releaseBlkPtr
  use Driver_interface, ONLY : Driver_abortFlash
  use gr_interface, ONLY : gr_GCAllocScratch,gr_GCReleaseScratch,&
       gr_GCTransferOneBlk

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent(IN) :: gridDataStruct
  logical, intent(IN) :: needSetup, releaseSetup
  integer, optional, intent(IN) :: blkCnt, indCnt
  integer, optional, dimension(:), intent(IN) :: blkList
  integer, optional, dimension(:), intent(IN) :: indList
  integer, optional, dimension(NDIM),intent(IN) :: gcCnt  
  logical, optional,intent(IN) :: putData
  
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  real, pointer, dimension(:,:,:,:) :: solnData
  integer :: blkID, ib

  integer :: locBlkCnt, locIndCnt
  real, pointer, dimension(:) :: locGC
  integer, pointer, dimension(:) :: locBlkList,locBlkoffsets,locIndList
  integer, dimension(1:NDIM) :: locGCCnt
  
  logical :: allPresent, getBlkInd, mode, locPutData
  
  if(releaseSetup)call gr_GCReleaseScratch(gridDataStruct)
  if(needSetup) then 
     allPresent = (present(blkCnt)).and.(present(blkList)).and.&
          (present(indCnt)).and.(present(indList)).and.present(gcCnt)
     if(.not.allPresent)call Driver_abortFlash("Grid_GCPutScratch: when needSetup is true all optional arguments must be present")
     call gr_GCAllocScratch(gridDataStruct,blkCnt,blkList,&
          indCnt,indList, gcCnt)
  end if

  locPutData=.true.
  if(present(putData)) locPutData=putData

  if(locPutData) then

     mode = .true.
  
     select case(gridDataStruct)
     case(CENTER)
        locBlkCnt = gr_GCCtrBlkCnt
        locIndCnt = gr_GCCtrIndCnt
        locGCCnt = gr_GCCtrGCCnt
        locBlkList => gr_GCCtrBlkList
        locIndList => gr_GCCtrIndList
        locBlkOffsets => gr_GCCtrBlkOffsets
        locGC => gr_GCCtr
     case(FACEX)
        locBlkCnt = gr_GCFxBlkCnt
        locIndCnt = gr_GCFxIndCnt
        locGCCnt = gr_GCFxGCCnt
        locBlkList => gr_GCFxBlkList
        locIndList => gr_GCFxIndList
        locBlkOffsets => gr_GCFxBlkOffsets
        locGC => gr_GCFx
     case(FACEY)
        locBlkCnt = gr_GCFyBlkCnt
        locIndCnt = gr_GCFyIndCnt
        locGCCnt = gr_GCFyGCCnt
        locBlkList => gr_GCFyBlkList
        locIndList => gr_GCFyIndList
        locBlkOffsets => gr_GCFyBlkOffsets
        locGC => gr_GCFy
     case(FACEZ)
        locBlkCnt = gr_GCFzBlkCnt
        locIndCnt = gr_GCFzIndCnt
        locGCCnt = gr_GCFzGCCnt
        locBlkList => gr_GCFzBlkList
        locIndList => gr_GCFzIndList
        locBlkOffsets => gr_GCFzBlkOffsets
        locGC => gr_GCFz
     end select
     
     do ib = 1,locBlkCnt
        blkID = locBlkList(ib)
        call Grid_getBlkIndexLimits(blkID,&
             blkLimits,blkLimitsGC,gridDataStruct)
        blkLimitsGC(LOW,1:NDIM)=blkLimits(LOW,1:NDIM)-locGCCnt(1:NDIM)
        blkLimitsGC(HIGH,1:NDIM)=blkLimits(HIGH,1:NDIM)+locGCCnt(1:NDIM)
        call Grid_getBlkPtr(blkID,solnData,gridDataStruct)
        call gr_GCTransferOneBlk(mode, &
          locIndCnt,locIndList,locBlkOffsets(ib),&
          blkLimits,blkLimitsGC,&
          locGC, solnData)
        call Grid_releaseBlkPtr(blkID,solnData,gridDataStruct)
     end do
  end if
end subroutine Grid_GCPutScratch
