!!****f* source/Grid/Grid_getBlkCenterCoords
!!
!! NAME
!!  Grid_getBlkCenterCoords
!!
!! SYNOPSIS
!!
!!  Grid_getBlkCenterCoords(integer(IN) :: blockId,
!!                         real(OUT)   :: blockCenter(MDIM))
!!                      
!!  
!! DESCRIPTION 
!!   Gets the coordinates of the center of the block identified by
!!   blockId.  Returns the coordinates in an array blockCenter
!!
!!
!! ARGUMENTS
!!  blockId - local id number of the block. for UG always 1
!!  blockCenter - returned array of size MDIM holding the blockCenter coords
!!
!! Example
!!   In 2 dimensions, if physical coordinates are ...
!!    
!!     ________________(0.5 1.0)
!!    |                |
!!    |                |
!!    |                |
!!    |                |
!!    |                |
!!    |                |
!!    |                |
!!    |_______________ |
!!  (-0.5, 0.0)
!!
!!  then the values returned in blockCenter are 
!!  blockCenter(IAXIS) = 0.0
!!  blockCenter(JAXIS) = 0.5
!!  blockCenter(KAXIS) = 0.0 since the dimension is not included  
!!
!!***


subroutine Grid_getBlkCenterCoords(blockId, blockCenter)

  implicit none
#include "constants.h"
  integer,intent(in) :: blockId
  real,dimension(MDIM),intent(out) :: blockCenter
  blockCenter=0.0
end subroutine Grid_getBlkCenterCoords
