!!****if* source/Grid/GridMain/paramesh/Grid_getCellCoords
!!
!! NAME
!!  Grid_getCellCoords
!!
!! SYNOPSIS
!!
!!  Grid_getCellCoords(integer(IN) :: axis,
!!                      integer(IN):: blockID, 
!!                      integer(IN):: edge, 
!!                      logical(IN):: guardcell, 
!!                      real(OUT)  :: coordinates(size),
!!                      integer(IN):: size)
!!  
!!  
!! DESCRIPTION
!!
!!    This subroutine is an accessor function that gets the coordinates of
!!    the cells in a given block.
!!    Coordinates are retrieved one axis at a time, 
!!    meaning you can get the i, j, _or_ k coordinates with one call.  
!!    If you want all the coordinates, all axes, you
!!    need to call Grid_getCellCoords 3 times, one for each axis.
!!    The code carries coordinates at cell centers as well as faces.
!!    It is possible to get coordinates for CENTER, only LEFT_EDGE,
!!    only RIGHT_EDGE or for all FACES along a dimension.
!!
!!
!!
!!
!! ARGUMENTS
!!            
!!   axis - specifies the integer index coordinates of the cells being retrieved.
!!          axis can have one of three different values, IAXIS, JAXIS or KAXIS 
!!          (defined in constants.h as 1,2 and 3)
!!
!!   blockID - integer block number
!!
!!   edge - integer value with one of four values, 
!!          LEFT_EDGE, RIGHT_EDGE, CENTER or FACES
!!          The edge argument specifies what side of the zone to get, 
!!          the CENTER point, the LEFT_EDGE  or the RIGHT_EDGE of the zone.
!!          FACES gets the left and right face of each cell, but since 
!!          two adjacent cells have a common face, there are only N+1
!!          unique values if N is the number of cells.
!!
!!   guardcell - logical input. If true coordinates for guardcells are returned
!!          along with the interior cells, if false, only the interior coordinates 
!!          are returned.
!!
!!          
!!   coordinates - The array holding the data returning the coordinate values
!!                 coordinates must be at least as big as "size" (see below)
!!           
!!   size - integer specifying the size of the coordinates array.
!!          if edge = CENTER/LEFT_EDGE/RIGHT_EDGE then
!!                If guardcell true then size =  interior cells + 2*guardcells
!!                otherwise size = number of interior cells
!!          If edge=FACES 
!!                If guardcell true then size =  interior cells + 2*guardcells+1
!!                otherwise size = number of interior cells+1
!!
!!               
!!  EXAMPLE 
!!
!!  1. Getting cell centered values
!!
!!   #include "constants.h"
!!   #include "Flash.h"
!!
!!      
!!      integer :: coordSize
!!      integer :: xCoord(coordSize) !sized to be number of coords returned
!!      
!!      do i=1, localNumBlocks
!!
!!          !get the index limits of the block
!!          call Grid_getBlkIndexLimits(i, blkLimits, blkLimitsGC)
!!
!!          !holds the number of cells returned in idir
!!          coordSize = blkLimitsGC(HIGH, IAXIS)
!!          call Grid_getCellCoords(IAXIS, i, CENTER, .true., xCoord, coordSize) 
!!
!!     end do    
!!
!!  2. Getting face values
!! 
!!   #include "constants.h"
!!   #include "Flash.h"
!!
!!      
!!      integer :: coordSize
!!      integer :: xCoord(coordSize) !sized to be number of coords returned
!!      
!!      do i=1, localNumBlocks
!!
!!          !get the index limits of the block
!!          call Grid_getBlkIndexLimits(i, blkLimits, blkLimitsGC)
!!
!!          !holds the number of cells returned in idir
!!          coordSize = blkLimitsGC(HIGH, IAXIS)+1
!!          call Grid_getCellCoords(IAXIS, i, FACES, .true., xCoord, coordSize) 
!!
!!     end do    
!!
!!
!!  NOTES
!!   variables that start with "gr_" are variables of Grid unit scope
!!   and are stored in the FORTRAN module Grid_data. Variables that are not
!!   starting with gr_ are local variables or arguments passed to the 
!!   routine.
!!
!!***

#ifdef DEBUG
#define DEBUG_GRID
#endif

subroutine Grid_getCellCoords(axis, blockID, edge, guardcell, coordinates, size)

  use Grid_data, ONLY : gr_oneBlock
  use Driver_interface, ONLY : Driver_abortFlash

#include "constants.h"
#include "Flash.h"

  implicit none

  integer, intent(in) :: axis,blockID, edge
  integer, intent(in) :: size
  logical, intent(in) :: guardcell
  real,intent(out), dimension(size) :: coordinates

  integer :: bOffset,eOffset,factor,calcSize,numGuard

  ! Do some error checking here
  

  if(((NDIM==1).and.(axis/=IAXIS)).or.((NDIM==2).and.(axis==KAXIS))) then
     bOffset = 0
     eOffset = 0
     numGuard = 0
  else
     if(guardcell) then
        bOffset = 0
        eOffset = 2*NGUARD
        numGuard = NGUARD
     else
        boffset = NGUARD
        eoffset = NGUARD
        numGuard = 0
     end if
  end if
  if(axis == IAXIS) then
     calcSize = NXB+2*numGuard
  else if(axis == JAXIS) then
     calcSize = NYB+2*K2D*numGuard
  else if(axis == KAXIS) then
     calcSize = NZB+2*K3D*numGuard
  end if

#ifdef DEBUG_GRID
  if (edge==FACES) calcSize = calcSize+1
  print*,' get coordinates', axis, blockID, edge, guardcell,size
  if((blockID<1).or.(blockID>MAXBLOCKS)) then
     call Driver_abortFlash("Grid_getCellCoords :invalid blockID ")
  end if
  if(.not.((edge==LEFT_EDGE).or.(edge==RIGHT_EDGE).or.&
       &(edge==CENTER))) then
     call Driver_abortFlash("Get Coords : invalid edge specification, must be LEFT_EDGE &
          RIGHT_EDGE, or CENTER")
  end if

!!!  This can be refined further to make it geometry specific

  if(.not.((axis==IAXIS).or.(axis==JAXIS).or.(axis==KAXIS))) then
     call Driver_abortFlash("Get Coords : invalid axis, must be IAXIS, JAXIS or KAXIS ")
  end if
  
  if(size < calcSize)then
     call Driver_abortFlash("Get Coords : size of output array too small")
  end if
#endif

  if(axis==IAXIS) then
     if(edge /= FACES) then
        coordinates(:) = &
             gr_oneBlock(BlockID)%firstAxisCoords&
             (edge,bOffset+1:eOffset+NXB)
     else
        coordinates(1:size-1)= gr_oneBlock(BlockID)%firstAxisCoords&
             (LEFT_EDGE,bOffset+1:eOffset+NXB)
        coordinates(size)=gr_oneBlock(BlockID)%firstAxisCoords&
             (RIGHT_EDGE,eOffset+NXB)
     end if
  elseif(axis==JAXIS) then
     if(edge /= FACES) then
        coordinates(:) = &
             gr_oneBlock(BlockID)%secondAxisCoords&
             (edge,bOffset+1:eOffset+NYB)
     else
        coordinates(1:size-1)= gr_oneBlock(BlockID)%secondAxisCoords&
             (LEFT_EDGE,bOffset+1:eOffset+NYB)
        coordinates(size)=gr_oneBlock(BlockID)%secondAxisCoords&
             (RIGHT_EDGE,eOffset+NYB)
     end if
  else
     if(edge /= FACES) then
        coordinates(:) = &
             gr_oneBlock(BlockID)%thirdAxisCoords(edge,bOffset+1:eOffset+NZB)
     else
        coordinates(1:size-1)= gr_oneBlock(BlockID)%thirdAxisCoords&
             (LEFT_EDGE,bOffset+1:eOffset+NZB)
        coordinates(size)=gr_oneBlock(BlockID)%thirdAxisCoords&
             (RIGHT_EDGE,eOffset+NZB)
     end if
  end if
  return
end subroutine Grid_getCellCoords





