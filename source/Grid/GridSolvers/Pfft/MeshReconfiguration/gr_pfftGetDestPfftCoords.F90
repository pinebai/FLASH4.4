!!****if* source/Grid/GridSolvers/Pfft/MeshReconfiguration/gr_pfftGetDestPfftCoords
!!
!! NAME
!!
!!  gr_pfftGetDestPfftCoords
!!
!! SYNOPSIS
!!
!!  gr_pfftGetDestPfftCoords(integer(1:MDIM), intent(IN) :: startPos, &
!!                           integer(1:MDIM), intent(OUT) :: pfftProcCoords)
!!
!! DESCRIPTION
!!
!!  Determine the coordinates of the PFFT processor that will
!!  be the first to receive data from the block at position, startPos.
!!  
!! ARGUMENTS
!!
!!  startPos: Starting position of the FLASH block in each dimension.
!!  pfftProcCoords:  The coordinates of the first PFFT processor that should
!!                  receive data from the FLASH block.
!! 
!!***

#include "Flash.h"
#include "constants.h"

subroutine gr_pfftGetDestPfftCoords(startPos, pfftProcCoords)

  use gr_pfftData, ONLY : pfft_procGrid
  use gr_pfftReconfigData, ONLY : pfft_procLookup
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none
  integer, dimension(1:MDIM), intent(IN) :: startPos
  integer, dimension(1:MDIM), intent(OUT) :: pfftProcCoords
  integer :: globalStartGridPoint, globalEndGridPoint, i, j

  pfftProcCoords = 0
  do i = 1, NDIM
     j = 0
     eachpfftProc: do        

        globalStartGridPoint = pfft_procLookup(i) % procInfo(j) % globalStartGridPoint
        globalEndGridPoint = pfft_procLookup(i) % procInfo(j) % globalEndGridPoint

        if ( (globalStartGridPoint <= startPos(i)).and.&
             (globalEndGridPoint >= startPos(i)) ) then 
           exit  !We found the relevant PFFT processor in this dimension.
        end if

        j = j + 1
        if (j >= pfft_procGrid(i)) then
           print *, "[gr_pfftGetDestPfftCoords]: No more processors available in this dimension"  
           call Driver_abortFlash("[gr_pfftGetDestPfftCoords]: Out of bounds in pfft_procLookup!")
        end if

     end do eachpfftProc
     pfftProcCoords(i) = j
  end do

end subroutine gr_pfftGetDestPfftCoords
