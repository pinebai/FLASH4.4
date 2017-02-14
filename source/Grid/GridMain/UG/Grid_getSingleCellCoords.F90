!!****if* source/Grid/GridMain/UG/Grid_getSingleCellCoords
!!
!! NAME
!!  Grid_getSingleCellCoords
!!
!! SYNOPSIS
!!
!!  Grid_getSingleCellCoords(integer(in)       ::  ind(MDIM), 
!!                           integer(in)       ::  blockid,
!!                           integer(in)       ::  edge,
!!                           integer(in)       ::  beginCount,
!!                           real(OUT)         ::  coords(MDIM))
!!  
!! DESCRIPTION 
!!
!!  Returns the coordinates of a single cell of a given block.
!!
!!
!! ARGUMENTS
!!
!!  ind - array holding the indices of the cell whose coordinates to return
!! 
!!  blockid - local block number
!!
!!  edge - indicates if user wants the left, center or right edge of the block
!!         options are LEFT_EDGE, RIGHT_EDGE, or CENTER.  These constants are
!!         defined in constants.h
!!
!!
!!  beginCount : tells the routine where to start index counting.  beginCount can
!!               be set to INTERIOR or EXTERIOR.  If INTERIOR is specified
!!               guardcell indices are not included and index 1 is the first interior cell. 
!!               If EXTERIOR is specified
!!               the first index, 1, is the left most guardcell.  See example
!!               below for more explanation.  (For most of the FLASH architecture code,
!!               we use EXTERIOR.  Some physics routines, however, find it helpful 
!!               only to work on the internal parts of the blocks (without
!!               guardcells) and wish to keep loop indicies  
!!               going from 1 to NXB without having to worry about finding 
!!               the correct offset for the number of guardcells.) 
!!               (INTERIOR and EXTERIOR are defined in constants.h)
!!
!!
!!    
!!  coords : returned coordinates of the specificed cell
!!
!!
!!
!!***

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif


subroutine Grid_getSingleCellCoords(ind, blockId,edge, beginCount,coords)

  use Grid_data,ONLY :gr_iCoords,gr_jCoords,gr_kCoords,gr_guard
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "constants.h"
#include "Flash.h"



  integer,dimension(MDIM), intent(in) :: ind
  integer, intent(in) :: blockId, edge
  integer, intent(in) :: beginCount
  real, dimension(MDIM), intent(out) :: coords

#ifdef DEBUG_GRID
  print*,' inside Grid_getSingleCellCoords', ind, blockId, edge, beginCount,coords
  if((blockid /= 1)) then
     print*,"Get Coords :invalid blockid "
     call Driver_abortFlash("Grid_getSingleCellCoords :invalid blockid ")
  end if
  if((edge/=LEFT_EDGE).and.(edge/=CENTER).and.(edge/=RIGHT_EDGE))&
       call Driver_abortFlash('Grid_getSingleCellCoords : invalid edge')

  print*, 'leaving the DEBUG_GRID statement'
#endif

  if(beginCount == EXTERIOR) then
     coords(IAXIS) = gr_iCoords(edge,ind(IAXIS),1)
     coords(JAXIS) = gr_jCoords(edge,ind(JAXIS),1)
     coords(KAXIS) = gr_kCoords(edge,ind(KAXIS),1)
  else if(beginCount == INTERIOR) then
     coords(IAXIS) = gr_iCoords(edge,ind(IAXIS)+gr_guard(IAXIS),1)
     coords(JAXIS) = gr_jCoords(edge,ind(JAXIS)+gr_guard(JAXIS),1)
     coords(KAXIS) = gr_kCoords(edge,ind(KAXIS)+gr_guard(KAXIS),1)
  else
     call Driver_abortFlash("Grid_getSingleCellCoords: beginCount has incorrect value")
  end if
  return
end subroutine Grid_getSingleCellCoords









