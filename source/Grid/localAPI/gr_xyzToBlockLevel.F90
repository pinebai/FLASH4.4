!!****if* source/Grid/localAPI/gr_xyzToBlockLevel
!!
!!  NAME
!!     gr_xyzToBlockLevel
!!
!!  SYNOPSIS
!!     call gr_xyzToBlockLevel(integer(in)  :: lev,
!!                              real  (in)  :: xyz(NDIM),
!!                             integer(OUT) :: ijk(NDIM))
!!
!!  DESCRIPTION
!!
!!    Maps spatial domain coordinates to block indexing coordinates.
!!    A block index coordinate is a 4-tuple of level number plus i,j,k indices.
!!    Level numbers begin at 1 for the coarsest level.  i,j,k indices begin at 0
!!    and go up to (inclusively) gr_nblock?*2^(level-1)-1 (with ? being X,Y,Z for
!!    i,j,k respectively).  Each level has twice the number of coords as its parent
!!    for each dimension.  This routine has no knowledge of the actual refinement
!!    pattern so it will gladly compute ijk's for blocks that dont exist or
!!    for levels that are beyond the maximum level of refinement.
!!
!!  ARGUMENTS
!!    lev: (in) 1-based level we want i,j,k block coordinates on
!!    xyz: (in) point in the domain we want mapped to a block (must be between gr_?min/gr_?max)
!!    ijk: (out) 0-based integer coordinate of block for the level requested
!!
!! NOTES
!!
!!  This is a stub implementation that does not do anything useful.
!!
!!  A more useful implementation is provided by PARAMESH Grid
!!  implementations.
!!***
subroutine gr_xyzToBlockLevel(lev, xyz, ijk)
#include "Flash.h"
  implicit none
  
  integer, intent(in) :: lev
  real, intent(in) :: xyz(NDIM)
  integer, intent(out) :: ijk(NDIM)
  
  ijk = 0

end subroutine
