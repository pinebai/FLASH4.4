!!****if* source/Grid/GridMain/UG/Grid_getBlkPhysicalSize
!!
!! NAME
!!  Grid_getBlkPhysicalSize
!!
!! SYNOPSIS
!!
!!  Grid_getBlkPhysicalSize(integer(IN)  :: blockId,
!!                          real(OUT) :: blockSize(MDIM))
!!  
!! DESCRIPTION 
!!  Gets the physical size of the block along each dimension by calculating
!!  the deltas of each direction times the number of zones.  It is important
!!  to calculate the size of a block consistently because rounding errors
!!  can occur in some cases.  Use this function to get the size of the 
!!  block rather than calculating it in a new way.  
!!
!! ARGUMENTS
!!  blockId -local block number
!!  blockSize - returned array of size MDIM holding the size of 
!!              each dimension of the block
!!
!! EXAMPLE
!!
!!    If there is 2-d block with 8 cells along IAXIS and 4 cells along
!!    JAXIS, and a single cell size (or delta) along each dimension is 0.1, 
!!    then the part of the physical domain occupied by the block is 
!!    0.8 x 0.4 in size, as in the block shown below
!!    
!!     - - - - - - - - (1.0,0.8)
!!    | | | | | | | | |
!!     - - - - - - - - 
!!    | | | | | | | | |
!!     - - - - - - - - 
!!    | | | | | | | | |
!!     - - - - - - - - 
!!    | | | | | | | | |
!!     - - - - - - - - 
!!   (0.2,0.4)
!!  
!!
!! NOTES
!!
!!  The physical block sizes returned by this routine are
!!  Calculated using deltas.
!!
!!  This interface is general enough to work even
!!  when the block sizes (in terms of number of cells) are not uniform.
!!  
!!
!!***


#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

subroutine Grid_getBlkPhysicalSize(blockId, blockSize)

  use Grid_data, ONLY :gr_iCoords,gr_jCoords,gr_kCoords
  use Grid_data, ONLY: gr_ilo,gr_ihi,gr_jlo,gr_jhi,gr_klo,gr_khi

  implicit none

#include "constants.h"
#include "Flash.h"
  integer,intent(in) :: blockId
  real,dimension(MDIM),intent(out) :: blockSize

  blockSize=0.0
  blockSize(IAXIS) = gr_iCoords(RIGHT_EDGE,gr_ihi,1) - &
                       gr_iCoords(LEFT_EDGE,gr_ilo,1)
  if(NDIM>1) blockSize(JAXIS) = gr_jCoords(RIGHT_EDGE,gr_jhi,1)-&
       gr_jCoords(LEFT_EDGE,gr_jlo,1)
  if(NDIM>2) blockSize(KAXIS) = gr_kCoords(RIGHT_EDGE,gr_khi,1)-&
       gr_kCoords(LEFT_EDGE,gr_klo,1)

  return
end subroutine Grid_getBlkPhysicalSize












