!!****if* source/Grid/GridSolvers/Pfft/MeshReconfiguration/PtToPt/gr_pfftFnArgCopyPencilToBuf
!!
!! NAME 
!!
!! gr_pfftFnArgCopyPencilToBuf
!!
!! SYNOPSIS
!!
!! gr_pfftFnArgCopyPencilToBuf(type(node), pointer  :: item)
!!
!! DESCRIPTION 
!!
!! Copies data from 1D pencil array into a 1D array within the passed node.
!!
!! ARGUMENTS
!!
!! item - pointer to a node.  The current implementation stores the nodes in 
!!        a linked list which includes ut_listMethods.includeF90 methods.
!!
!! NOTES
!!
!! References pfft_outArray module level array.
!! A return value of 0 indicates successful copy.
!!
!!***

integer function gr_pfftFnArgCopyPencilToBuf(item)
#include "constants.h"
  use ut_conversionInterface, ONLY : ut_convertToMemoryOffset
  use gr_pfftNodeObject, ONLY : node
  use gr_pfftData, ONLY : pfft_inLen
  use gr_pfftReconfigData, ONLY : pfft_pfftBuf
  implicit none
  type(node), pointer  :: item
  integer :: index, i, j, k
  integer :: memoryOffset, oneDimensionalCoord
  integer, dimension(MDIM) :: elementCoord, arrayLBound, copyStart, copyEnd

  arrayLBound(:) = 1  !Required for ut_convertToMemoryOffset subroutine.
  index = 1
  copyStart(1:MDIM) = item % metadata % pfftStartPos
  copyEnd(1:MDIM) = item % metadata % pfftEndPos
  
  !pfft_inLen is the equal amount of space assigned to each PFFT processor.
  do k = copyStart(KAXIS), copyEnd(KAXIS); elementCoord(KAXIS) = k
     do j = copyStart(JAXIS), copyEnd(JAXIS); elementCoord(JAXIS) = j
        do i = copyStart(IAXIS), copyEnd(IAXIS); elementCoord(IAXIS) = i
           call ut_convertToMemoryOffset(MDIM, elementCoord, arrayLBound, &
                pfft_inLen, memoryOffset)
           oneDimensionalCoord = memoryOffset + 1  !Fortran arrays 1-based.
           item % buf(index) = pfft_pfftBuf(oneDimensionalCoord)
           index = index + 1
        end do
     end do
  end do
  gr_pfftFnArgCopyPencilToBuf = 0
end function gr_pfftFnArgCopyPencilToBuf
