!!****f* source/Grid/Grid_getSingleCellVol
!!
!! NAME
!!  Grid_getSingleCellVol
!!
!! SYNOPSIS
!!
!!  Grid_getSingleCellVol(integer(IN) :: blockID,
!!                        integer(IN) :: beginCount,
!!                        integer(IN) :: point(MDIM), 
!!                        real(OUT)   :: cellVolume)
!!  
!! DESCRIPTION 
!!  
!!  Gets cell volumes for a single cell in a given block. 
!!  
!! ARGUMENTS 
!!
!!  blockID - integer local blockid
!!
!!  beginCount : tells the routine where to start index counting.  beginCount can
!!               be set to INTERIOR or EXTERIOR.  If INTERIOR is specified,
!!               guardcell indices are not included and index 1 is the first interior cell. 
!!               If EXTERIOR is specified,
!!               the first index, 1, is the leftmost guardcell.  See example
!!               below for more explanation.  (For most of the FLASH architecture code,
!!               we use EXTERIOR.  Some physics routines, however, find it helpful 
!!               only to work on the internal parts of the blocks (without
!!               guardcells) and wish to keep loop indicies  
!!               going from 1 to NXB without having to worry about finding 
!!               the correct offset for the number of guardcells.) 
!!               (INTERIOR and EXTERIOR are defined in constants.h)
!! 
!!  point(MDIM):
!!           specifies the point to return
!!   
!!           point(1) = i
!!           point(2) = j
!!           point(3) = k
!!
!!           If a problem is only 2d, point(3) is ignored.  For 1d problems
!!           point(2) and point(3) are ignored.
!!
!!
!!  cellVolume - real value containing the cell volume
!!
!!
!! NOTES
!! 
!!  Current implementations of this interface assume that all cells in a 
!!  dimension of a block have the same grid spacing. The grid spacings used
!!  are the ones returned by Grid_getDeltas.
!!
!! SEE ALSO
!!
!!  Grid_getDeltas
!!  Grid_getSingleCellCoords
!!
!!***

subroutine Grid_getSingleCellVol(blockID, beginCount, point, cellvolume)

implicit none
#include "constants.h"

  integer, intent(in) :: blockID, beginCount
  integer, intent(in) :: point(MDIM)
  real, intent(out)   :: cellvolume

  cellvolume = 0.0
  return
end subroutine Grid_getSingleCellVol
