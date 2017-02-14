!!****if* source/Grid/GridSolvers/Pfft/MeshReconfiguration/Collective/gr_pfftHandleKaxisFragments
!!
!! NAME 
!!
!!   gr_pfftHandleKaxisFragments
!!
!! SYNOPSIS
!!
!!   gr_pfftHandleKaxisFragments(integer (IN) :: direction)
!!
!! DESCRIPTION 
!!
!! Subroutine handles the KAXIS communication.  It accepts a direction
!! argument which specifies if we are going to PFFT grid or 
!! returning from PFFT grid.
!!
!!  1. Disassembles the JAXIS block and places the KAXIS block 
!!     fragments into the appropriate send buffer slot.  The data in 
!!     the send buffer will end up on the correct KAXIS PFFT processor
!!     after the MPI_Alltoall.
!!       - required by Grid_pfftMapToInput. (FORWARDS)
!!
!!  2. Reassembles the KAXIS block fragments into complete 
!!     JAXIS blocks.  We will place the JAXIS block data into the
!!     appropriate send buffer slot.  The data in 
!!     the send buffer will end up on the correct FLASH processor
!!     after the MPI_Alltoall.
!!       - required by Grid_pfftMapFromOutput. (BACKWARDS)
!!
!!  The MPI_Alltoall is called immediately after this subroutine.
!!
!! ARGUMENTS
!!
!!   direction - indicates whether we are going to PFFT grid, or 
!!               returning from PFFT grid.
!!
!! NOTES 
!! 
!! There is the option to view the communication buffers by compiling
!! with DEBUG_COMM_BUFFERS.  This is designed to be used with the 
!! DEBUG_MAPPING mode in Grid_solvePoisson in which we communicate 
!! simple integer data.  Hence, the subroutine PrintCommunicationBuffer
!! casts the data to integer first for ease of viewing.
!!
!!***

subroutine gr_pfftHandleKaxisFragments(direction)

#include "Flash.h"
#include "constants.h"
#include "Pfft.h"

  use gr_pfftReconfigData, ONLY : pfft_maxProcs, pfft_sendKMap, pfft_recvJMap, &
       pfft_sendKMap, pfft_sendBuf, pfft_recvBuf
  use gr_pfftData, ONLY : pfft_procGrid, pfft_me, &
       pfft_ndim, pfft_commWithTopology
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, Grid_getBlkCornerID, &
       Grid_getBlkIndexLimits
  use Driver_interface, ONLY : Driver_abortFlash
  use gr_pfftInterface, ONLY : gr_pfftPrintCommBuffers
  use Logfile_interface, ONLY : Logfile_open, Logfile_close
  use Grid_data, ONLY : gr_meshMe

  implicit none

#include "Flash_mpi.h"

  integer, intent(IN) :: direction
  integer, allocatable, dimension(:) :: kMapPtr,sendPtr
  integer :: j, nblks, kb, ke
  integer :: recvInd, sendInd, recvMapInd, sendMapInd
  integer :: bsize, bufSize, blk, error
  integer :: destIndex, totalBlockSize, blockSizePacked, destProc, ierr, logUnit
  integer, dimension(1:MDIM) :: pfftProcCoords
  logical, parameter :: logUnitLocal = .true.

  if ( (direction /= TO_PFFT) .and. (direction /= FROM_PFFT) ) then
     call Driver_abortFlash("[gr_pfftHandleKaxisFragments]: Direction not recognised!")
  end if


#ifdef DEBUG_COMM_BUFFERS
  if (direction == FROM_PFFT) then
     !Data has just been received FROM the KAXIS PFFT processor(s).
     call Logfile_open(logUnit,logUnitLocal)
     write(logUnit,*) ""
     write(logUnit,*) "[Data from PFFT grid to PFFT grid (K to J Movement)] (cast to int)"
     call gr_pfftPrintCommBuffers(pfft_recvBuf, pfft_sendKMap, logUnit)
     call Logfile_close(logUnitLocal)
  end if
#endif


  !We must now extract the relevant data from pfft_recvBuf, and then 
  !send this to the correct KAXIS processor.
  allocate(kMapPtr(pfft_maxProcs), sendPtr(pfft_maxProcs), STAT=error)
  if (error /= 0) then
     call Driver_abortFlash("[gr_pfftHandleKaxisFragments]: Severe error. Memory cannot be allocated!")
  end if

  kMapPtr = PFFT_MAPEND
  sendPtr = 0

  do j = 1,pfft_maxProcs
     recvInd = 0
     recvMapInd = PFFT_MAPEND

     !This is the number of blocks we have just received after the first MPI_alltoall.
     nblks = pfft_recvJMap(PFFT_NUMBLKS,j)

     do blk = 1,nblks

        !Extract information about the single block.
        !kb and ke contain the original block coordinates in FLASH grid.
        kb = pfft_recvJMap(recvMapInd+PFFT_SPOSK,j)
        ke = pfft_recvJMap(recvMapInd+PFFT_EPOSK,j)

        totalBlockSize = pfft_recvJMap(recvMapInd+PFFT_BSIZE,j)
        blockSizePacked = 0
        bsize = pfft_recvJMap(recvMapInd+PFFT_BSIZE,j) / (ke-kb+1)
        pfftProcCoords(1:MDIM) = pfft_me(1:MDIM)  !i.e. me.

        !We cannot know the number of K block fragments per 
        !J block because of the way the data structure is defined.  Instead, we must
        !iterate through the K blocks until all the J block data is packed.
        do while (blockSizePacked < totalBlockSize)
           
           if (pfftProcCoords(KAXIS) >= pfft_procGrid(KAXIS)) then
              print *, "[gr_pfftHandleKaxisFragments]: Out of bounds"
              call Driver_abortFlash("[gr_pfftHandleKaxisFragments]: PFFT processor does not exist!")
           end if

           call MPI_Cart_rank(pfft_commWithTopology, pfftProcCoords(1), destProc, ierr)
           destIndex = destProc + 1

           sendMapInd = kMapPtr(destIndex)
           sendInd = sendPtr(destIndex)

           !kb and ke now contain the destination block coordinates in pfft grid.
           kb = pfft_sendKMap(sendMapInd+PFFT_SPOSK,destIndex)
           ke = pfft_sendKMap(sendMapInd+PFFT_EPOSK,destIndex)
           bufSize = bsize * (ke-kb+1)    


           if (direction == TO_PFFT) then
              !We are breaking the block in pfft_recvBuf into smaller
              !fragments that we will send to the appropriate KAXIS
              !processor.
              pfft_sendBuf(sendInd+1:sendInd+bufSize, destIndex) = &
                   pfft_recvBuf(recvInd+1:recvInd+bufSize, j)

           else if (direction == FROM_PFFT) then              
              !We are joining the fragments from the KAXIS PFFT 
              !processor into a single block which we will send
              !to the appropriate JAXIS PFFT processor.
              pfft_sendBuf(recvInd+1:recvInd+bufSize, j) = &
                   pfft_recvBuf(sendInd+1:sendInd+bufSize, destIndex)

           end if


           !Update pointer to the next section of kmap.
           kMapPtr(destIndex) = kMapPtr(destIndex) + PFFT_MAPEND

           !Update the cumulative count of data being sent to each proc.
           sendPtr(destIndex) = sendPtr(destIndex) + bufSize
           recvInd = recvInd + bufSize


           blockSizePacked = blockSizePacked + bufSize
           !Assert that we have not tried to send more data than we received!
           if (blockSizePacked > totalBlockSize) then
              print *, "[gr_pfftHandleKaxisFragments]: Block size violation."
              call Driver_abortFlash("[gr_pfftHandleKaxisFragments]: Block size violation.")
           end if


           pfftProcCoords(KAXIS) = pfftProcCoords(KAXIS) + 1     !Skip to next KAXIS processor.

        end do

        recvMapInd = recvMapInd + PFFT_MAPEND
     end do
  end do

  deallocate(kMapPtr, sendPtr, STAT=error)
  if (error /= 0) then
     call Driver_abortFlash("[Grid_pfftMapToInput]: Severe error. Memory cannot be deallocated!")
  end if


#ifdef DEBUG_COMM_BUFFERS
  if (direction == TO_PFFT) then
     !Data will be communicated TO the KAXIS PFFT processor(s).
     call Logfile_open(logUnit,logUnitLocal)
     write(logUnit,*) ""
     write(logUnit,*) "[Data from PFFT grid to PFFT grid (J to K Movement)] (cast to int)"
     call gr_pfftPrintCommBuffers(pfft_sendBuf, pfft_sendKMap, logUnit)
     call Logfile_close(logUnitLocal)
  end if
#endif

  return
end subroutine gr_pfftHandleKaxisFragments
