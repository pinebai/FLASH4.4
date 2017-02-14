!!****if* source/Grid/GridSolvers/Pfft/MeshReconfiguration/Collective/gr_pfftGenMap
!!
!! NAME
!!
!!  gr_pfftGenMap
!!
!! SYNOPSIS
!!  
!!  gr_pfftGenMap()
!!
!! DESCRIPTION
!!
!!  This subroutine generates a map from AMR or UG to Pfft
!!  The same routine will work for both meshes, because it does not
!!  assume that one block moves without fragmentation
!!  The algorithm works as follows. 
!!
!!  1. On processor (j1,k1) identify all data that needs to move to processor 
!!  (j2,:) for all values of j2. So all data meant for 
!!  processors (j2,1), (j2,2) ... (j2,n) is moved to (j2,k1) through an 
!!  alltoall along j axis
!!
!!  2. On processor (j2,k1), separate data meant for (j2,k2) for all values
!!  of k2 and send it to appropriate destination through and alltoall along
!!  k axis.
!!
!! ARGUMENTS
!!
!! NOTES
!!
!! This routine does not send the actual data, but uses communication
!! to generate the information about the communication pattern to be expected
!! when actual data movement happens.
!!
!!***

subroutine gr_pfftGenMap()

#include "Flash.h"
#include "constants.h"
#include "Pfft.h"

  use Logfile_interface, ONLY : Logfile_stamp
  use gr_pfftData, ONLY : pfft_procGrid, pfft_ndim, pfft_globalLen, &
       pfft_comm, pfft_inLen, pfft_myPE, pfft_regionBndBox,&
       pfft_inRegion
  use gr_pfftReconfigData, ONLY : pfft_maxProcData, pfft_maxProcs, &
       pfft_sendJMap, pfft_recvJMap, pfft_sendKMap, pfft_recvKMap, &
       pfft_procLookup, pfft_fragmentPtr, pfft_sendBuf, pfft_recvBuf
  use gr_pfftinterface, ONLY : gr_pfftGenMapHelper, &
       gr_pfftCopyToSendMap, gr_pfftGridPointTable
  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkCornerID,&
       Grid_getBlkIndexLimits
  use Driver_interface, ONLY : Driver_abortFlash
#ifdef FLASH_GRID_PARAMESH
  use Grid_data, ONLY : gr_oneRefLev
#endif


  implicit none

#include "Flash_mpi.h"

  integer,dimension(IAXIS:KAXIS) :: stride,cornerID,startPos,endPos,&
       blkSize,destCoords
  integer,dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  integer, dimension(MAXBLOCKS) :: blkList
  integer, dimension(MDIM) :: maxOverlap
  integer :: i,j,k,destProc,destIndex,blockID
  integer :: jEnd,kEnd,kb,ke,bSize,ierr
  integer :: blkCount, error, mapStart, mapSize
  integer :: localMapSize, maxJProcData, maxKProcData, axis, pfftProc
  integer :: totalBlockFragmentCount, numProcsToReceiveBlocks, numOverlap
  integer :: firstFlashOverlap, pfftProcStart, pfftProcEnd, flashStartPosition

  if(pfft_ndim == 1) then
     return
  end if

  !Modifies pfft_procLookup.
  call gr_pfftGridPointTable(pfft_inLen)


#if defined(FLASH_GRID_UG)
  call Grid_getListOfBlocks(LEAF,blkList,blkCount)

#elif defined(FLASH_GRID_PARAMESH)

  if(pfft_inRegion) then
     if(gr_oneRefLev==NONEXISTENT) then
        call Grid_getListOfBlocks(INREGION,blkList,blkCount,region_bndBox=pfft_regionBndBox)        
     else
        call Grid_getListOfBlocks(INREGION,blkList,blkCount,gr_oneRefLev,pfft_regionBndBox)        
     end if
  else
     call Grid_getListOfBlocks(REFINEMENT,blkList,blkCount,gr_oneRefLev)
  end if

#endif


  !The first thing we must do is calculate the maximum number of 
  !block fragments that will be communicated in the MPI_Alltoall.
  !--------------------------------------------------------------
  blockID = blkList(1)  !We just need 1 block as all blocks are at same level.
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  blkSize(1:MDIM) = blkLimits(HIGH,1:MDIM) - blkLimits(LOW,1:MDIM) + 1

  maxOverlap(:) = 1
  do axis = 1, pfft_ndim

     !blockEnd is the actual end destination of the block, and pfftProcEnd is the
     !end destination of the PFFT grid processor.
     do pfftProc = 0, pfft_procGrid(axis)-1

        numOverlap = 0
        pfftProcStart = pfft_procLookup(axis) % procInfo(pfftProc) % globalStartGridPoint
        pfftProcEnd = pfft_procLookup(axis) % procInfo(pfftProc) % globalEndGridPoint

        !Now I know the coordinates of the PFFT processor, I need to find 
        !the number of FLASH blocks that occupy the same space.
        firstFlashOverlap = (pfftProcStart / blkSize(axis))
        flashStartPosition = 1 + (firstFlashOverlap * blkSize(axis))

        do while (flashStartPosition <= pfftProcEnd)
           numOverlap = numOverlap + 1           
           flashStartPosition = flashStartPosition + blkSize(axis)
        end do

        maxOverlap(axis) = max(maxOverlap(axis), numOverlap)
     end do
  end do

  totalBlockFragmentCount = product(maxOverlap)
  !--------------------------------------------------------------
  !We multiply the total number of block fragments by the amount of 
  !data used to describe the block fragment.  We also add 
  !a single block fragment multiplied by the data used to 
  !desribe that block fragment.  This extra section is 
  !is used to store aggregate quantities.
  localMapSize = PFFT_MAPEND * (totalBlockFragmentCount + 1)


  call MPI_ALLREDUCE(localMapSize,mapSize,1,FLASH_INTEGER,MPI_MAX,&
       pfft_comm(IAXIS),ierr)
  !--------------------------------------------------------------


  !Allocate the fragment pointer array which is used to keep track
  !of where we currently are in the "map".  Also allocate the JAXIS 
  !"maps" which are used to describe the data to be sent and received
  !along the JAXIS.
  !--------------------------------------------------------------
  pfft_maxProcs = product(pfft_procGrid(1:NDIM))

  allocate(pfft_fragmentPtr(pfft_maxProcs), STAT=error)
  if (error /= 0) then
     call Driver_abortFlash("Severe error: Memory cannot be allocated!")
  end if

  allocate(pfft_sendJMap(mapSize,pfft_maxProcs), &
       pfft_recvJMap(mapSize,pfft_maxProcs), &
       STAT=error)
  if (error /= 0) then
     call Driver_abortFlash("Severe error: Memory cannot be allocated!")
  end if

  pfft_sendJMap = -5; pfft_recvJMap = -5 !Temp initialisation: Easier to spot errors.
  pfft_fragmentPtr(:) = PFFT_MAPEND  !Skip past the aggregate data block.
  !--------------------------------------------------------------


  !Determine which processor(s) each block will be sent to.  If the block 
  !belongs to several processors in PFFT grid, then we will send a 
  !fragment of the block to each processor.  We first consider the 
  !JAXIS movement.
  !--------------------------------------------------------------
  do i = 1, blkCount
     blockID = blkList(i)
     call Grid_getBlkCornerID(blockID,cornerID,stride)
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     blkSize(1:MDIM) = blkLimits(HIGH,1:MDIM) - blkLimits(LOW,1:MDIM) + 1

     !Starting position of the block if this was assumed to be the 
     !UG at least refined level
     startPos(1:MDIM) = 1
     startPos(1:NDIM) = ceiling(real(cornerID(1:NDIM)) / real(stride(1:NDIM)))
     endPos(1:MDIM) = 1
     endPos(1:NDIM) = startPos(1:NDIM) + blkSize(1:NDIM) - 1
     bSize = product(blkSize(1:NDIM)) / blkSize(JAXIS)

     call gr_pfftCopyToSendMap(blockID, startPos, endPos, bsize, JAXIS)

  end do
  !Calculate the aggregate quantities for the JAXIS
  call gr_pfftGenMapHelper(JAXIS, pfft_fragmentPtr, maxJProcData)
  pfft_maxProcData = maxJProcData
  !--------------------------------------------------------------


  !MPI treats pfft_sendJMap as a contiguous buffer divided into sections of mapSize.
  !In Fortran abstraction, the second dimension is of size pfft_maxProcs. 
  !Therefore the 3rd "mapSize" block is destined for processor 3.
  call MPI_ALLTOALL(pfft_sendJMap(1,1),mapSize,FLASH_INTEGER,&
       pfft_recvJMap(1,1),mapSize,FLASH_INTEGER,pfft_comm(IAXIS),ierr)


  !call PrintSendMapInfo(JAXIS,pfft_maxProcs)

  !At this point all the information about what to expect from 
  !other processors has been obtained. Now we need to separate the
  !Z direction communication from the overall picture.
  !The 1st data transfer involved the movement from FLASH grid to PFFT grid 
  !along the JAXIS.  Now we handle the 2nd data transfer phase along the KAXIS.
  !--------------------------------------------------------------
  if(pfft_ndim > 2) then

     allocate(pfft_sendKMap(mapSize,pfft_maxProcs), &
          pfft_recvKMap(mapSize,pfft_maxProcs), &         
          STAT=error)
     if (error /= 0) then
        call Driver_abortFlash("Severe error. Memory cannot be allocated!")
     end if

     pfft_fragmentPtr(:) = PFFT_MAPEND  !Reset pfft_fragmentPtr.

     !Loop over all received metadata from each processor.
     do i = 1, pfft_maxProcs
        j = PFFT_MAPEND

        !Loop over each received block from the 1st data transfer phase.
        do k = 1,pfft_recvJMap(PFFT_NUMBLKS,i)

           blockID = pfft_recvJMap(j+PFFT_BLKID,i)
           startPos = pfft_recvJMap(j+PFFT_SPOSI:j+PFFT_SPOSK,i)
           endPos = pfft_recvJMap(j+PFFT_EPOSI:j+PFFT_EPOSK,i)
           bSize = pfft_recvJMap(j+PFFT_BSIZE,i) / (endPos(KAXIS)-startPos(KAXIS)+1)

           call gr_pfftCopyToSendMap(blockID, startPos, endPos, bsize, KAXIS)

           j = j + PFFT_MAPEND  !Skip to next received block.

        end do  !Loop over each received block from a single processor.
     end do  !Loop over each processor.

     !Calculate the aggregate quantities for the KAXIS
     call gr_pfftGenMapHelper(KAXIS, pfft_fragmentPtr, maxKProcData)
     pfft_maxProcData = max(pfft_maxProcData, maxKProcData)

     !call PrintSendMapInfo(KAXIS,pfft_maxProcs)
  
     call MPI_ALLTOALL(pfft_sendKMap(1,1), mapSize, FLASH_INTEGER,&
          pfft_recvKMap(1,1), mapSize, FLASH_INTEGER, pfft_comm(IAXIS), ierr)

  end if
  !--------------------------------------------------------------

  deallocate(pfft_fragmentPtr, STAT=error)
  if (error /= 0) then
     call Driver_abortFlash("Severe error. Memory cannot be deallocated!")
  end if


  !Finally allocate the send/receive buffers.
  if (pfft_maxProcData > 0) then
     !There is the same allocation on each processor in this implementaion.
     !This allocation is for the persistent communication buffers.
     !*2 for send & recv buffers and *8 for 8-byte reals.
     call Logfile_stamp( pfft_maxProcData*pfft_maxProcs*2*8, &
          "[gr_pfftGenMap] Actual comm. buffer allocation (bytes)")

     allocate(pfft_sendBuf(pfft_maxProcData,pfft_maxProcs), &
          pfft_recvBuf(pfft_maxProcData,pfft_maxProcs), STAT=error)
     if (error /= 0) then
        call Driver_abortFlash &
             ("[gr_pfftGenMap]: Severe error. Memory cannot be allocated!")
     end if
  end if

  return
end subroutine gr_pfftGenMap


!For debugging only.
subroutine PrintSendMapInfo(axis, numProcs)

#include "Flash.h"
#include "constants.h"
#include "Pfft.h"

  use gr_pfftData, ONLY : pfft_myPE
  use gr_pfftReconfigData, ONLY : pfft_sendJMap, pfft_sendKMap
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none
  integer, intent(IN) :: axis, numProcs
  integer :: nblks, i
  integer, dimension(:,:), pointer :: sendMap

  !We always query the receive map regardless of the direction.
  if (axis == JAXIS) then
     sendMap => pfft_sendJMap
  else if (axis == KAXIS) then
     sendMap => pfft_sendKMap
  else
     call Driver_abortFlash("[PrintSendMapInfo]: Axis must be JAXIS or KAXIS!")
  end if
  !-----------------------------------------------------------------------

  do i = 1, numProcs
     nblks = sendMap(PFFT_NUMBLKS,i)
     if (nblks /= 0) then
        print *, "FORWARD... Processor", pfft_myPE, "axis:", axis, "sending", nblks, "to processor", i-1
     end if
  end do

end subroutine PrintSendMapInfo
