!!****if* source/Grid/GridSolvers/Pfft/MeshReconfiguration/Collective/gr_pfftPrintCommBuffers
!!
!! NAME
!!
!!  gr_pfftPrintCommBuffers
!!
!! SYNOPSIS
!!
!!  gr_pfftPrintCommBuffers (real, dimension(:,:), intent(IN) :: buffer, &
!!                           integer, dimension(:,:), intent(IN) :: map, &
!!                           integer, intent(IN) :: logUnit)
!!
!! DESCRIPTION
!!
!!  This subroutine is to help with debugging.  It simply prints the
!!  communication buffer on each processor.  Results are cast to int because
!!  it is intended to be used with the DEBUG mode in the customised 
!!  Grid_solvePoisson.
!!
!! ARGUMENTS
!!
!!  buffer:  Specifies the communication buffer to print out.
!!  map:  Specifies the map which contains the layout of the communication buffer.
!!  logUnit:  The file handle of an already open file.
!!
!!***
subroutine gr_pfftPrintCommBuffers(buffer, map, logUnit)

#include "Flash.h"
#include "constants.h"
#include "Pfft.h"

  use gr_pfftReconfigData, ONLY : pfft_maxProcs 
  use gr_pfftData, ONLY : pfft_myPE

  implicit none
  real, dimension(:,:), intent(IN) :: buffer
  integer, dimension(:,:), intent(IN) :: map
  integer, intent(IN) :: logUnit
  integer,dimension(IAXIS:KAXIS) :: startCoords, endCoords
  integer :: i, startElement, endElement, offset, numBlocks, blk

  !We want to iterate over each slot in the communication buffer.
  do i = 1, pfft_maxProcs

     startElement = 1
     offset = PFFT_MAPEND
     numBlocks = map(PFFT_NUMBLKS,i)  !can contain 0 blocks
     write(logUnit,*) numBlocks, "blocks assigned to processor", i-1    

     do blk = 1, numBlocks
        
        startCoords(1:MDIM) = map(offset+PFFT_SPOSI:offset+PFFT_SPOSK,i)
        endCoords(1:MDIM) = map(offset+PFFT_EPOSI:offset+PFFT_EPOSK,i)

        endElement = startElement + ( &
             (endCoords(KAXIS) - startCoords(KAXIS) + 1) * &
             (endCoords(JAXIS) - startCoords(JAXIS) + 1) * &
             (endCoords(IAXIS) - startCoords(IAXIS) + 1) - 1 )

        write(logUnit,*) "Send block", blk, "PFFT coords i:", startCoords(IAXIS), endCoords(IAXIS), &
             ", j:", startCoords(JAXIS), endCoords(JAXIS), ", k:", &
             startCoords(KAXIS), endCoords(KAXIS)
        write(logUnit,*) "Send block", blk, "data:", int(buffer(startElement:endElement, i))

        startElement = endElement + 1        
        offset = offset + PFFT_MAPEND

     end do
  end do

end subroutine gr_pfftPrintCommBuffers
