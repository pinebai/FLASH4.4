!!****if* source/Grid/GridMain/UG/Grid_getDeltas
!!
!! NAME
!!  Grid_getDeltas
!!
!! SYNOPSIS
!!
!!  call Grid_getDeltas(integer(IN):: blockId, 
!!                      real(OUT)::del(MDIM))
!!  
!! DESCRIPTION 
!!  
!!  Gets the dx/dy/dz for a given blockId on the Grid
!!  dx is the size of one zone in the x direction of a block
!!
!!  
!! ARGUMENTS 
!!
!!  blockId - local block number
!!  del - array of size MDIM returned holding the dx, dy, and dz values
!!
!!  
!!***
#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

subroutine Grid_getDeltas(blockId,del)
  use Grid_data, ONLY : gr_delta
  implicit none
#include "constants.h"

  integer, intent(IN) :: blockId
  real, dimension(MDIM),intent(OUT) :: del
  
  del = gr_delta
  return
end subroutine Grid_getDeltas

