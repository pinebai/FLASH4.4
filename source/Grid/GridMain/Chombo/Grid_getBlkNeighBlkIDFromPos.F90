!!****if* source/Grid/GridMain/Chombo/Grid_getBlkNeighBlkIDFromPos
!!
!! NAME
!!
!!  gr_findNeghID
!!
!! SYNOPSIS
!!
!!  call Grid_getBlkNeighBlkIDFromPos(
!!                            real(IN)    :: blockID, 
!!                            real(IN)    :: pos(MDIM), 
!!                            integer(IN) :: neghDir(MDIM),
!!                            integer(OUT) :: ansBlockID,
!!                            integer(OUT) :: ansProcID)
!!
!!  call Grid_getBlkIDFromPos(real(IN)    :: blockID, 
!!                            real(IN)    :: pos(MDIM), 
!!                            integer(IN) :: neghDir(MDIM),
!!                            integer(OUT) :: ansBlockID,
!!                            integer(OUT) :: ansProcID)
!!
!! DESCRIPTION
!!
!!   Given the physical coordinates of a point outside the current
!!   block but in its neighborhood, and the direction in which the
!!   neighboring block lies, his routine finds the processor number and blockID
!!   within that processor for the neighboring block that contains 
!!   the point.
!!
!! ARGUMENTS
!!
!!   blockID : ID of block in current processor
!!
!!   pos :     coordinates of the point of interest
!!
!!   neghDir : the location of the neighbor with respect to the
!!             current block, in other words specification on which
!!             face/edge/point is common between the current block and
!!             neighbor of interest. For example
!!             negh(1:MDIM)=LEFT_EDGE indicates that the lowest
!!             left hand corner of the current block is the same as
!!             the highest right end corner of the neighbor. Similarly
!!             negh(IAXIS)=LEFT_EDGE, negh(JAXIS:KAXIS) =
!!             CENTER implies that the left face of current block is
!!             common with the right face of the neighbor.
!!             This array is in the form returned by Grid_outsideBoundBox.
!!
!!   ansBlockID : identity of the neighbor, the first number 
!!                is the blocknumber within the processor
!!   ansProcID :  identity of the neighbor, the second number is the processor
!!                number where the neighbor is located
!!
!! NOTES
!!
!!   Currently only implemented for PARAMESH Grid implementations.
!!
!!   The specific subroutine Grid_getBlkNeighBlkIDFromPos is also available
!!   under the generic name Grid_getBlkIDFromPos.
!!
!! SEE ALSO
!!   Grid_outsideBoundBox
!!***

subroutine Grid_getBlkNeighBlkIDFromPos(blockID,pos,neghDir,ansBlockID,ansProcID)

  use Driver_interface,ONLY: Driver_abortFlash
  
  implicit none

#include "constants.h"
  integer,intent(IN) :: blockID
  real,dimension(MDIM),intent(IN) :: pos
  integer,dimension(MDIM),intent(IN) :: neghDir
  integer,intent(OUT) :: ansBlockID, ansProcID

  call Driver_abortFlash("Grid_getBlkNeighBlkIDFromPos only implemented for PARAMESH Grids!")
  ansProcID=0
  ansBlockID=0

end subroutine Grid_getBlkNeighBlkIDFromPos
