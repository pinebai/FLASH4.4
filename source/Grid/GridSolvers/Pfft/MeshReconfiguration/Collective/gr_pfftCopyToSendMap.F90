!!****if* source/Grid/GridSolvers/Pfft/MeshReconfiguration/Collective/gr_pfftCopyToSendMap
!!
!! NAME
!!
!!  gr_pfftCopyToSendMap
!!
!! SYNOPSIS
!!
!!  gr_pfftCopyToSendMap(integer(IN)   :: blockID,
!!                       integer(IN) :: blockStartPos(MDIM),
!!                       integer(IN) :: blockEndPos(MDIM),
!!                       integer(IN)   :: bsize,
!!                       integer(IN)   :: axis)
!!
!! DESCRIPTION
!!
!!  This routine is used to distribute a block (given by blockID) of global 
!!  position (given by blockStartPos and blockEndPos) to the processors 
!!  in PFFT grid.  
!!
!!  This routine will be called once for the JAXIS distribution and once 
!!  for the KAXIS distribution for a particular block.  The behaviour is 
!!  slightly different for JAXIS distribtion and KAXIS distribution:
!!
!!  * When we are working with the JAXIS we must locate the first PFFT 
!!    processor which overlaps with the FLASH block.  We then distribute the block 
!!    to this processor and (where necessary) subsequent processors along the same
!!    axis.
!!  * When we are working with the KAXIS we have just received blocks 
!!    from the JAXIS distribution.  As such, the first PFFT processor along 
!!    the KAXIS is simply this processor's ID.  Where necessary, the block
!!    will be distributed to subsequent processors along the KAXIS.
!!
!!
!! ARGUMENTS
!!
!!  blockID       - ID of the block to be distributed.
!!  blockStartPos - Array containing start position of the block along each axis.
!!  blockEndPos   - Array containing end position of the block along each axis.
!!  bsize         - The size of the block (or sometimes block fragment if KAXIS
!!                  distribution) that must be distributed to PFFT processors.
!!  axis          - Specifies whether this is JAXIS or KAXIS distribution.
!!
!!***

subroutine gr_pfftCopyToSendMap(blockID, blockStartPos, blockEndPos, bsize, axis)

#include "Flash.h"
#include "constants.h"
#include "Pfft.h"

  use gr_pfftData, ONLY : pfft_procGrid, pfft_ndim, pfft_globalLen, pfft_me, &
       pfft_commWithTopology
  use gr_pfftReconfigData, ONLY : pfft_sendJMap, pfft_sendKMap, &
       pfft_maxProcs, pfft_fragmentPtr, pfft_procLookup
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getBlkCornerID, Grid_getBlkIndexLimits
  use gr_pfftinterface, ONLY : gr_pfftGetDestPfftCoords

  implicit none
  integer, intent(IN) :: blockID
  integer, dimension(1:MDIM), intent(IN) :: blockStartPos, blockEndPos
  integer, intent(IN) :: bsize, axis

  integer, dimension(1:MDIM) :: startPos, endPos
  integer, dimension(:,:), pointer :: sendMap
  integer :: destProc, destIndex, mapStart, pfftProcEnd, blockEnd, ierr
  integer, dimension(1:MDIM) :: pfftProcCoords


  startPos(1:MDIM) = blockStartPos(1:MDIM)
  endPos(1:MDIM) = blockEndPos(1:MDIM)
  pfftProcCoords(1:MDIM) = 0

  if (axis == JAXIS) then
     sendMap => pfft_sendJMap

     !Search for the coordinates of the PFFT processor which will 
     !be the first destination for the current block.
     call gr_pfftGetDestPfftCoords(startPos, pfftProcCoords)

  else if (axis == KAXIS) then
     sendMap => pfft_sendKMap
     pfftProcCoords(1:MDIM) = pfft_me(1:MDIM)

  else
     call Driver_abortFlash("[gr_pfftCopyToSendMap]: Axis must be JAXIS or KAXIS!")
  end if

  blockEnd = endPos(axis)


  !Loop over each processor along this axis
  eachProcAlongAxis: do

     if (pfftProcCoords(axis) >= pfft_procGrid(axis)) then
        print *, "[gr_pfftCopyToSendMap]: Out of bounds"
        call Driver_abortFlash("[gr_pfftCopyToSendMap]: PFFT processor does not exist!")
     end if    

     call MPI_Cart_rank(pfft_commWithTopology, pfftProcCoords(1), destProc, ierr)
     destIndex = destProc + 1


     !mapStart is an integer which indicates the next available space in our 
     !stored data structure.
     mapStart = pfft_fragmentPtr(destIndex)

     !blockEnd is the actual end destination of the block, and pfftProcEnd is the
     !end destination of the PFFT grid processor.
     pfftProcEnd = pfft_procLookup(axis) % procInfo(pfftProcCoords(axis)) % globalEndGridPoint


     if (pfftProcEnd < blockEnd) then
        !The block end point exists beyond this PFFT processor.
        endPos(axis) = pfftProcEnd
     else
        endPos(axis) = blockEnd
     end if


     !Update the map and fragmentPtr with the block fragment's details.
     sendMap(mapStart+PFFT_BLKID,destIndex) = blockID
     sendMap(mapStart+PFFT_SPOSI:mapStart+PFFT_SPOSK,destIndex) = startPos
     sendMap(mapStart+PFFT_EPOSI:mapStart+PFFT_EPOSK,destIndex) = endPos
     sendMap(mapStart+PFFT_BSIZE,destIndex) = bSize * (endPos(axis)-startPos(axis)+1)
     pfft_fragmentPtr(destIndex) = pfft_fragmentPtr(destIndex) + PFFT_MAPEND

     !Break when the entire block has been distributed amongst this axis processors.
     if (endPos(axis) == blockEnd) then
        exit
     else if (pfftProcEnd == pfft_globalLen(axis)) then
        !We may have some blocks with zero size which should not be touched.
        !We would expect the block to be distributed before hitting zero sized blocks.
        print *, "[gr_pfftCopyToSendMap]: Block not yet distributed to processors"
        call Driver_abortFlash("[gr_pfftCopyToSendMap]: Block distribution error!")
     end if


     !If we are here, the end criterion was not met so send a fragment to 
     !the next processor.
     startPos(axis) = endPos(axis) + 1
     pfftProcCoords(axis) = pfftProcCoords(axis) + 1

  end do eachProcAlongAxis

  nullify(sendMap)

end subroutine gr_pfftCopyToSendMap
