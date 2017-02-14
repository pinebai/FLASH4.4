!!****if* source/Grid/GridMain/UG/Grid_getCellCoords
!!
!! NAME
!!  Grid_getCellCoords
!!
!! SYNOPSIS
!!
!!  Grid_getCellCoords(integer(IN)  :: axis,
!!                      integer(IN) :: blockId, 
!!                      integer(IN) :: edge, 
!!                      logical(IN) :: guardcell, 
!!                      real(OUT)   :: coordinates(size),
!!                      integer(IN) :: size)
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
!!   blockId - integer block number
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
!!   and are stored in the fortran module Grid_data. Variables are not
!!   starting with gr_ are local varibles or arguments passed to the 
!!   routine.
!!
!!***


#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif


subroutine Grid_getCellCoords(axis, blockId, edge, guardcell,coordinates, size)
  use Grid_data, ONLY : gr_guard,gr_iCoords,gr_jCoords,gr_kCoords,&
       gr_ilo,gr_ihi,gr_jlo,gr_jhi,gr_klo,gr_khi
  use Driver_interface, ONLY : Driver_abortFlash

#include "constants.h"
#include "Flash.h"

  implicit none
  integer, intent(in) :: axis
  integer, intent(in) :: blockId, edge
  integer, intent(in) :: size
  logical, intent(in) :: guardcell
  real,intent(out), dimension(size) :: coordinates

  integer :: bOffset,eOffset,numGuard,calcSize


  bOffset=0
  eOffset=0
  if(axis <= NDIM) then
     if(guardcell) then
        bOffset = 0
        eOffset = 2*gr_guard(axis)
        numGuard = gr_guard(axis)
     else
        boffset = gr_guard(axis)
        eoffset = gr_guard(axis)
        numGuard = 0
     end if
  end if
#ifdef DEBUG_GRID
  print*,' inside Grid_getCellCoords', axis, blockId, edge, numGuard,size
  if((blockid /= 1)) then
     print*,"Get Coords :invalid blockid "
     call Driver_abortFlash("Get Coords :invalid blockid ")
  end if
  if(.not.((edge==LEFT_EDGE).or.(edge==RIGHT_EDGE).or.&
       &(edge==CENTER))) then
     print*,"Get Coords : invalid edge specification"
     call Driver_abortFlash("Get Coords : invalid edge specification")
  end if
  
!!!  This can be refined further to make it geometry specific
  if(.not.((axis==IAXIS).or.(axis==JAXIS).or.(axis==KAXIS))) then
     print*,"Get Coords : invalid axis "
     call Driver_abortFlash("Get Coords : invalid axis ")
  end if
  

  if(axis==IAXIS)calcSize=gr_ihi-gr_ilo+1+2*numGuard
  if(axis==JAXIS)calcSize=gr_jhi-gr_jlo+1+2*numGuard*K2D
  if(axis==KAXIS)calcSize=gr_khi-gr_klo+1+2*numGuard*K3D
  if (edge==FACES) calcSize = calcSize+1
  if(size < calcSize)then
     print*,"Get Coords : size of output array too small",size,calcSize
     call Driver_abortFlash("Get Coords : size of output array too small")
  end if

!  print*, 'leaving the DEBUG_GRID statement'
#endif

  
  if(axis==IAXIS) then
     if(edge /=FACES) then
        coordinates(1:size) = gr_iCoords(edge,&
             bOffset+1:eOffset+gr_ihi-gr_ilo+1,1)
     else
        coordinates(1:size-1)=gr_iCoords(LEFT_EDGE,&
             bOffset+1:eOffset+gr_ihi-gr_ilo+1,1)
        coordinates(size)=gr_iCoords(RIGHT_EDGE,eoffset+gr_ihi-gr_ilo+1,1)
     end if
  elseif(axis==JAXIS) then
     if(edge /=FACES) then
        coordinates(1:size) = gr_jCoords(edge,&
             bOffset+1:eOffset+gr_jhi-gr_jlo+1,1)
     else
        coordinates(1:size-1)=gr_jCoords(LEFT_EDGE,&
             bOffset+1:eOffset+gr_jhi-gr_jlo+1,1)
        coordinates(size)=gr_jCoords(RIGHT_EDGE,eoffset+gr_jhi-gr_jlo+1,1)
     end if
  else
     if(edge /=FACES) then
        coordinates(1:size) = gr_kCoords(edge,&
             bOffset+1:eOffset+gr_khi-gr_klo+1,1)
     else
        coordinates(1:size-1)=gr_kCoords(LEFT_EDGE,&
             bOffset+1:eOffset+gr_khi-gr_klo+1,1)
        coordinates(size)=gr_kCoords(RIGHT_EDGE,eoffset+gr_khi-gr_klo+1,1)
     end if
  end if
  return
end subroutine Grid_getCellCoords
