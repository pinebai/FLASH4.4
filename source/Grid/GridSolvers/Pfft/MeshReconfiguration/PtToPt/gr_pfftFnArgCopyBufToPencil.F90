!!****if* source/Grid/GridSolvers/Pfft/MeshReconfiguration/PtToPt/gr_pfftFnArgCopyBufToPencil
!!
!! NAME 
!!
!! gr_pfftFnArgCopyBufToPencil
!!
!! SYNOPSIS
!!
!! gr_pfftFnArgCopyBufToPencil(type(node), pointer  :: item)
!!
!! DESCRIPTION 
!!
!! Copies data from a 1D array within the passed node to the 
!! final 1D pencil array.
!!
!! ARGUMENTS
!!
!! item - pointer to a node.  The current implementation stores the nodes in 
!!        a linked list which includes ut_listMethods.includeF90 methods.
!!
!! SIDE EFFECTS
!!
!! Inserts values into pfft_inArray.
!!
!! NOTES
!!
!! A return value of 0 indicates successful copy.
!!
!!***

integer function gr_pfftFnArgCopyBufToPencil(item)
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
           pfft_pfftBuf(oneDimensionalCoord) = item % buf(index)
           index = index + 1
        end do
     end do
  end do
  gr_pfftFnArgCopyBufToPencil = 0
end function gr_pfftFnArgCopyBufToPencil
