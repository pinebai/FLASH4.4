!!****if* source/Grid/GridMain/paramesh/Grid_getSingleCellCoords
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

  use Grid_data, ONLY : gr_oneBlock
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "constants.h"
#include "Flash.h"



  integer,dimension(MDIM), intent(in) :: ind
  integer, intent(in) :: blockId, edge
  integer, intent(in) :: beginCount
  real, dimension(MDIM), intent(out) :: coords


#ifdef DEBUG_GRID
  print*,' inside Grid_getSingleCellCoords', ind, blockId, edge, beginCount
#ifdef FLASH_GRID_UG
  if((blockid /= 1)) then
     print*,"Grid_getSingleCellCoords :invalid blockid "
     call Driver_abortFlash("Grid_getSingleCellCoords :invalid blockid ")
  end if
#endif
  if((edge/=LEFT_EDGE).and.(edge/=CENTER).and.(edge/=RIGHT_EDGE))&
       call Driver_abortFlash('Grid_getSingleCellCoods : invalid edge')
  if(beginCount == EXTERIOR) then
     if((ind(IAXIS)<GRID_ILO_GC).or.(ind(IAXIS)>GRID_IHI_GC))&
          call Driver_abortFlash('GetSingleCellCoords : I index out of blkLimits')
     if((ind(JAXIS)<GRID_JLO_GC).or.(ind(JAXIS)>GRID_JHI_GC))&
          call Driver_abortFlash('GetSingleCellCoords : J index out of blkLimits')
     if((ind(KAXIS)<GRID_KLO_GC).or.(ind(KAXIS)>GRID_KHI_GC))&
          call Driver_abortFlash('GetSingleCellCoords : K index out of blkLimits')
  else if(beginCount == INTERIOR) then
     if((ind(IAXIS)<GRID_ILO).or.(ind(IAXIS)>GRID_IHI))&
          call Driver_abortFlash('GetSingleCellCoords : I index out of blkLimits')
     if((ind(JAXIS)<GRID_JLO).or.(ind(JAXIS)>GRID_JHI))&
          call Driver_abortFlash('GetSingleCellCoords : J index out of blkLimits')
     if((ind(KAXIS)<GRID_KLO).or.(ind(KAXIS)>GRID_KHI))&
          call Driver_abortFlash('GetSingleCellCoords : K index out of blkLimits')
  else
     call Driver_abortFlash("Grid_getSingleCellCoords, incorrect value for beginCount")
  end if

!  print*, 'leaving the DEBUG_GRID statement'
#endif


  if(beginCount == EXTERIOR) then
     coords(IAXIS) = gr_oneBlock(BlockId)%firstAxisCoords(edge,ind(IAXIS))
     coords(JAXIS) = gr_oneBlock(BlockId)%secondAxisCoords(edge,ind(JAXIS))
#if NDIM==3
     coords(KAXIS) = gr_oneBlock(BlockId)%thirdAxisCoords(edge,ind(KAXIS))
#endif
  else if(beginCount == INTERIOR) then
     coords(IAXIS)=gr_oneBlock(BlockId)%firstAxisCoords(edge,ind(IAXIS)+NGUARD)
     coords(JAXIS)=gr_oneBlock(BlockId)%secondAxisCoords(edge,ind(JAXIS)+NGUARD)
#if NDIM==3
     coords(KAXIS)=gr_oneBlock(BlockId)%thirdAxisCoords(edge,ind(KAXIS)+NGUARD)
#endif
  endif


  return
end subroutine Grid_getSingleCellCoords









