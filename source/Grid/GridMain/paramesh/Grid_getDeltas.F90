!!****if* source/Grid/GridMain/paramesh/Grid_getDeltas
!!
!! NAME
!!  Grid_getDeltas
!!
!! SYNOPSIS
!!
!!  Grid_getDeltas(integer(IN) :: blockId,
!!                 real(OUT)   :: del(MDIM))
!!  
!! DESCRIPTION 
!!  
!!  Gets the grid spacing dx/dy/dz for a given blockId on the Grid.
!!  dx is the size of one cell in the x direction of a block.
!!  
!!  
!! ARGUMENTS 
!!
!!  blockId - local block number
!!  del - array of size MDIM returned holding the dx, dy, and dz values
!!
!!***

subroutine Grid_getDeltas(blockId,del)

  use tree, ONLY : lrefine
  use Grid_data, ONLY : gr_delta

  implicit none

#include "constants.h"

  integer, intent(IN)   :: blockId
  real, dimension(MDIM), intent(out) :: del
  
  del = gr_delta(1:MDIM,lrefine(blockId))
  return
end subroutine Grid_getDeltas
