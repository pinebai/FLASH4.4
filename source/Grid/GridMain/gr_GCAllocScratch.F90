!!****if* source/Grid/GridMain/gr_GCAllocScratch
!!
!! NAME
!!  gr_GCAllocScratch
!!
!! SYNOPSIS
!!
!!  gr_GCAllocScratch(integer(IN) :: gridDataStruct,
!!                    integer(IN) :: blkCnt,
!!                    integer(IN) :: blkList(blkCnt),
!!                    integer(IN) :: indCnt,
!!                    integer(IN) :: indList(indCnt),
!!                    integer(IN) :: gcCnt(NDIM))
!!  
!!  
!! DESCRIPTION
!!  
!!  This routine sets up space for storing 
!!  guardcells of mesh data structures such as cell centered 
!!  or face centered variables of a given list of blocks. It is possible
!!  to specify the count of guardcells to be saved along each dimension,
!!  and also the specific variable indices into the corresponding mesh
!!  data structures. 
!!
!!  This routine needs to be called once for a given collection of blocks,
!!  data structure, and indices. The corresponding get and put functions 
!!  can be called several times as long as there is no change in blocks 
!!  and indices for a data structure. The accompanying 
!!  gr_GCReleaseScratch function releases the space allocated in 
!!  this function.
!!
!!  If there is any change in the mesh, list of blocks, or the number
!!  or order of indices, the gr_releaseGuardScratchSpace function must be 
!!  called, followed by a new call to gr_GCAllocScratch
!!  with appropriately changed input arguments.
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
!!  blkCnt     : count of blocks whose guard cells are to be saved
!!  blkList    : list of blocks whose guard cells are to be saved
!!  indCnt     : count of indices in the data structure to be saved
!!  indList    : list of indices in the data structure to be saved
!!  gcCnt      : count of guardcells along each dimension that need to be saved                  
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
!!   gr_GCReleaseScratch
!!   Grid_GCPutScratch and Grid_GCTransferOneBlk
!!
!!***

#include "constants.h"
#include "Flash.h"

subroutine gr_GCAllocScratch(gridDataStruct,blkCnt,blkList,&
     indCnt,indList, gcCnt)
  
  use gr_GCScratchData

  use Grid_interface, ONLY : Grid_getBlkIndexLimits

  implicit none

  integer, intent(IN) :: gridDataStruct, blkCnt, indCnt
  integer, dimension(blkCnt), intent(IN) :: blkList
  integer, dimension(indCnt), intent(IN) :: indList
  integer, dimension(NDIM), intent(IN) :: gcCnt

  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  integer :: arraysize, blkId, blkSize, interiorSize
  integer :: i,j,k, blkIndCnt
  integer, allocatable, dimension(:) :: blkOffsets

  allocate(blkOffsets(blkCnt))


  arraysize=0
  do i=1,blkCnt
     blkId=blkList(i)
     call Grid_getBlkIndexLimits(blkID,blkLimits,blkLimitsGC,gridDataStruct)
     blkSize=1; interiorSize=1
     do j = 1,NDIM
        blkIndCnt = blkLimits(HIGH,j)-blkLimits(LOW,j)+1
        interiorSize=interiorSize*blkIndCnt
        blkSize=blkSize*(blkIndCnt + 2*gcCnt(j))
     end do
     blkOffsets(i)=arraysize
     arraysize=arraysize+indCnt*(blksize-interiorSize)
  end do
  
  select case(gridDataStruct)
  case(CENTER)
     allocate(gr_GCCtrIndlist(indCnt))
     allocate(gr_GCCtrBlkList(blkCnt))
     allocate(gr_GCCtrBlkOffsets(blkCnt))
     allocate(gr_GCCtr(arraysize))
     gr_GCCtrIndCnt=indCnt
     gr_GCCtrblkCnt=blkCnt
     gr_GCCtrBlkOffsets(1:blkCnt)=blkOffsets(1:blkCnt)
     gr_GCCtrGCCnt(1:NDIM) = gcCnt(1:NDIM)
  case(FACEX)
     allocate(gr_GCFxIndlist(indCnt))
     allocate(gr_GCFxBlkList(blkCnt))
     allocate(gr_GCFxBlkOffsets(blkCnt))
     allocate(gr_GCFx(arraysize))
     gr_GCFxIndCnt=indCnt
     gr_GCFxblkCnt=blkCnt
     gr_GCFxBlkOffsets(1:blkCnt)=blkOffsets(1:blkCnt)
     gr_GCFxGCCnt(1:NDIM) = gcCnt(1:NDIM)
  case(FACEY)
     allocate(gr_GCFyIndlist(indCnt))
     allocate(gr_GCFyBlkList(blkCnt))
     allocate(gr_GCFyBlkOffsets(blkCnt))
     allocate(gr_GCFy(arraysize))
     gr_GCFyIndCnt=indCnt
     gr_GCFyblkCnt=blkCnt
     gr_GCFyBlkOffsets(1:blkCnt)=blkOffsets(1:blkCnt)
     gr_GCFyGCCnt(1:NDIM) = gcCnt(1:NDIM)
  case(FACEZ)
     allocate(gr_GCFzIndlist(indCnt))
     allocate(gr_GCFzBlkList(blkCnt))
     allocate(gr_GCFzBlkOffsets(blkCnt))
     allocate(gr_GCFz(arraysize))
     gr_GCFzIndCnt=indCnt
     gr_GCFzblkCnt=blkCnt
     gr_GCFzBlkOffsets(1:blkCnt)=blkOffsets(1:blkCnt)
     gr_GCFzGCCnt(1:NDIM) = gcCnt(1:NDIM)
  end select
  
  deallocate(blkOffsets)
  
end subroutine gr_GCAllocScratch
