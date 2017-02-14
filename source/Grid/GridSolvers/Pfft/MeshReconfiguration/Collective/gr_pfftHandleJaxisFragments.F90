!!****if* source/Grid/GridSolvers/Pfft/MeshReconfiguration/Collective/gr_pfftHandleJaxisFragments
!!
!! NAME 
!!
!!   gr_pfftHandleJaxisFragments
!!
!! SYNOPSIS
!!
!!   gr_pfftHandleJaxisFragments(integer (IN) :: direction, &
!!                               integer (IN) :: gridVar)
!!
!! DESCRIPTION 
!!
!! Subroutine handles the JAXIS communication.  It accepts a direction
!! argument which specifies if we are going to PFFT grid or 
!! returning from PFFT grid.
!!  
!! 1. Copies data from the grid into the appropriate slot in 
!!    pfft_sendBuf.  After the MPI_Alltoall, the data will end up 
!!    on the correct JAXIS PFFT processor. (FORWARDS)
!! 2. Copies data from pfft_recvBuf into the appropriate location
!!    in the grid. (BACKWARDS)
!!
!! ARGUMENTS
!!
!!   direction - indicates whether we are going to PFFT grid, or 
!!               returning from PFFT grid.
!!   gridVar   - variable on the mesh on which pfft is to be applies
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

!!REORDER(4): solnData

subroutine gr_pfftHandleJaxisFragments(direction, gridVar)

#include "Flash.h"
#include "constants.h"
#include "Pfft.h"

  use gr_pfftReconfigData, ONLY : pfft_maxProcs, pfft_sendJMap, pfft_sendBuf, pfft_recvBuf
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, Grid_getBlkCornerID, &
       Grid_getBlkIndexLimits
  use Driver_interface, ONLY : Driver_abortFlash
  use gr_pfftInterface, ONLY : gr_pfftPrintCommBuffers
  use Logfile_interface, ONLY : Logfile_open, Logfile_close
  use Grid_data, ONLY : gr_meshMe

  implicit none

#include "Flash_mpi.h"

  integer, intent(IN) :: direction, gridVar

  integer,dimension(IAXIS:KAXIS) :: stride, cornerID, startCoords, endCoords, &
       blkSize, blockLocation, guard
  integer,dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  real, dimension(:,:,:,:), pointer :: solnData
  integer :: blockID, i, ii, jj, kk, offset, blk, index, logUnit
  logical, parameter :: logUnitLocal = .true.

  if ( (direction /= TO_PFFT) .and. (direction /= FROM_PFFT) ) then
     call Driver_abortFlash("[gr_pfftHandleJaxisFragments]: Direction not recognised!")
  end if


#ifdef DEBUG_COMM_BUFFERS
  if (direction == FROM_PFFT) then
     !Data has just been received FROM the JAXIS PFFT processor(s).
     call Logfile_open(logUnit,logUnitLocal)
     !If we only do data movement and no FFT, the data is the same as TO_PFFT.
     write(logUnit,*) ""
     write(logUnit,*) "[Data from PFFT grid (J-Movement) to FLASH grid] (cast to int)"
     call gr_pfftPrintCommBuffers(pfft_recvBuf, pfft_sendJMap, logUnit)
     call Logfile_close(logUnitLocal)
  end if
#endif

  !pfft_sendJMap describes which data is sent from FLASH processors
  !to JAXIS PFFT processors.
  do i = 1,pfft_maxProcs
     index = 1
     offset = PFFT_MAPEND

     do blk = 1,pfft_sendJMap(PFFT_NUMBLKS,i)  !can contain 0 blocks

        !! From the saved map find the block number
        blockID = pfft_sendJMap(offset+PFFT_BLKID,i)

        call Grid_getBlkCornerID(blockID,cornerID,stride)
        call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
        guard(1:MDIM) = blkLimits(LOW,1:MDIM) - blkLimitsGC(LOW,1:MDIM)
        blkSize(1:MDIM) = blkLimits(HIGH,1:MDIM) - blkLimits(LOW,1:MDIM) + 1

        call Grid_getBlkPtr(blockID,solnData,CENTER)

        !We need to convert global grid point identifiers to local grid point identifiers.
        blockLocation(1:MDIM) = 0
        blockLocation(1:NDIM) = (cornerID(1:NDIM)-1) / stride(1:NDIM)

        startCoords(1:MDIM) = pfft_sendJMap(offset+PFFT_SPOSI:offset+PFFT_SPOSK,i) - &
             blockLocation(1:MDIM) + guard(1:MDIM)
        endCoords(1:MDIM) = pfft_sendJMap(offset+PFFT_EPOSI:offset+PFFT_EPOSK,i) - &
             blockLocation(1:MDIM) + guard(1:MDIM)


        !Loops are contained inside "if" statements for efficiency.
        if (direction == TO_PFFT) then
           do kk = startCoords(KAXIS), endCoords(KAXIS)
              do jj = startCoords(JAXIS), endCoords(JAXIS)
                 do ii = startCoords(IAXIS), endCoords(IAXIS)
                    pfft_sendBuf(index,i) = solnData(gridVar,ii,jj,kk)
                    index = index + 1
                 end do
              end do
           end do
        else if (direction == FROM_PFFT) then
           do kk = startCoords(KAXIS), endCoords(KAXIS)
              do jj = startCoords(JAXIS), endCoords(JAXIS)
                 do ii = startCoords(IAXIS), endCoords(IAXIS)
                    solnData(gridVar,ii,jj,kk) = pfft_recvBuf(index,i)
                    index = index + 1
                 end do
              end do
           end do
        end if


        call Grid_releaseBlkPtr(blockID,solnData,CENTER)
        offset = offset + PFFT_MAPEND
     end do

  end do

#ifdef DEBUG_COMM_BUFFERS
  if (direction == TO_PFFT) then
     !Data will be communicated TO the JAXIS PFFT processor(s).
     call Logfile_open(logUnit,logUnitLocal)
     write(logUnit,*) ""
     write(logUnit,*) "[Data from FLASH grid to PFFT grid (J-Movement)] (cast to int)"
     call gr_pfftPrintCommBuffers(pfft_sendBuf, pfft_sendJMap, logUnit)
     call Logfile_close(logUnitLocal)
  end if
#endif

  return
end subroutine gr_pfftHandleJaxisFragments
