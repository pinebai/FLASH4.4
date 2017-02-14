!!****if* source/Grid/GridMain/UG/Grid_getBlkCenterCoords
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
!!  
!! DESCRIPTION 
!!  Gets the coordinates of the center of the block identified by
!!  blockId.  Returns the coordinates in an array blockCenter  
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
!!
!!***

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif


subroutine Grid_getBlkCenterCoords(blockId,blockCenter)
  use Grid_interface, ONLY : Grid_getBlkBoundBox

  implicit none

#include "constants.h"
#include "Flash.h"
  integer,intent(in) :: blockId
  real,dimension(MDIM),intent(out) :: blockCenter
  integer :: i

  real,dimension(2,MDIM) :: bndBox
  call Grid_getBlkBoundBox(blockId, bndBox)

  blockCenter = 0.0

  blockCenter(1) = bndBox(1,1) + (bndBox(2,1) - bndBox(1,1))/2
  if(NDIM>1) blockCenter(2) = bndBox(1,2) + (bndBox(2,2) - bndBox(1,2))/2
  if(NDIM>2) blockCenter(3) = bndBox(1,3) + (bndBox(2,3) - bndBox(1,3))/2


  return
end subroutine Grid_getBlkCenterCoords














