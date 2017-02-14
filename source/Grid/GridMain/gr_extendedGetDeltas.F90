!!****if* source/Grid/GridMain/gr_extendedGetDeltas
!!
!! NAME
!!  gr_extendedGetDeltas
!!
!! SYNOPSIS
!!
!!  gr_extendedGetDeltas(integer(IN):: blockID, 
!!                       integer(IN):: pe, 
!!                       real(OUT)  :: del(MDIM))
!!  
!!  
!! DESCRIPTION
!!
!!    This subroutine is an accessor function that gets the grid spacing
!!    of the cells in a given block. This is a generalized variant of
!!    Grid_getDeltas, for use within the Grid unit only:
!!    it can be called for blocks local to the executing processor or,
!!    at least in the Paramesh3 implementation of the Grid unit, for
!!    remote blocks for which cached information is currently locally
!!    available.
!!
!!    The block about which information is requested is identified by
!!    a pair (blockID, pe). All other arguments are as for Grid_getDeltas
!!    and are passed to Grid_getDeltas if it is called by this routine.
!!    This is a generic implementation which always fails when (pe .NE. myPE).
!!
!!
!!
!!
!! ARGUMENTS
!!            
!!   blockID - integer block number
!!
!!   pe      - processor where block (or a cached copy of block info) resides
!!
!!  del      - array of size MDIM returned holding the dx, dy, and dz values
!!               
!!
!!
!!  NOTES
!!   variables that start with "gr_" are variables of Grid unit scope
!!   and are stored in the fortran module Grid_data. Variables are not
!!   starting with gr_ are local variables or arguments passed to the 
!!   routine.
!!
!!  SEE ALSO
!!
!!   Grid_getDeltas
!!
!!***

#ifdef DEBUG
#define DEBUG_GRID
#endif

subroutine gr_extendedGetDeltas(blockID, pe, del)

  use Grid_data, ONLY : gr_meshMe
  use Grid_interface, ONLY : Grid_getDeltas
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none
#include "constants.h"

  integer, intent(in) :: blockID,pe
  real, dimension(MDIM),intent(OUT) :: del


  if (pe.EQ.gr_meshMe) then
     call Grid_getDeltas(blockID, del)
  else
     call Driver_abortFlash('Calling gr_extendedGetDeltas for uncached'// & 
             &      ' remote blocks is not supported in this Grid implementation.')
  endif

  return
end subroutine gr_extendedGetDeltas





