!!****if* source/Grid/GridMain/Chombo/Grid_getCellCoords
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

subroutine Grid_getCellCoords(axis, blockId, edge, guardcell,coordinates, size)
  use chombo_f_c_interface, ONLY : ch_get_cell_coords
  implicit none
  integer, intent(in) :: axis
  integer, intent(in) :: blockId, edge
  integer, intent(in) :: size
  logical, intent(in) :: guardcell
  real,intent(out), dimension(size) :: coordinates
  integer :: gCell

  if (guardcell) then
     gCell = 1
  else
     gCell = 0
  end if
  call ch_get_cell_coords(blockId, axis, edge, size, gCell, coordinates)
end subroutine Grid_getCellCoords
