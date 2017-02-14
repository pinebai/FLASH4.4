!!****if* source/Grid/GridSolvers/Pfft/MeshReconfiguration/PtToPt/gr_pfftFnArgCopyGridToBuf
!!
!! NAME 
!!
!! gr_pfftFnArgCopyGridToBuf
!!
!! SYNOPSIS
!!
!! gr_pfftFnArgCopyGridToBuf(type(node), pointer  :: item)
!!
!! DESCRIPTION 
!!
!! Copies data from the grid into a 1D array within the passed node.
!!
!! ARGUMENTS
!!
!! item - pointer to a node.  The current implementation stores the nodes in 
!!        a linked list which includes ut_listMethods.includeF90 methods.
!!
!! NOTES
!!
!! References pfft_srcGridVar module level variable.
!! A return value of 0 indicates successful copy.
!!
!!***

!!REORDER(4): solnData

integer function gr_pfftFnArgCopyGridToBuf(item)
#include "constants.h"
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr
  use gr_pfftNodeObject, ONLY : node
  use gr_pfftReconfigData, ONLY : pfft_srcGridVar
  implicit none
  type(node), pointer  :: item
  integer, dimension(MDIM) :: startCoords, endCoords
  real, dimension(:,:,:,:), pointer :: solnData
  integer :: blockID, gridVar, ii, jj, kk, index

  gridVar = pfft_srcGridVar  !Module level variable.
  index = 1
  blockID = item % metadata % flashBlockID
  startCoords(1:MDIM) = item % metadata % flashStartPos(1:MDIM)
  endCoords(1:MDIM) = item % metadata % flashEndPos(1:MDIM)

  call Grid_getBlkPtr(blockID,solnData,CENTER)
  do kk = startCoords(KAXIS), endCoords(KAXIS)
     do jj = startCoords(JAXIS), endCoords(JAXIS)
        do ii = startCoords(IAXIS), endCoords(IAXIS)
           item % buf(index) = solnData(gridVar,ii,jj,kk)
           index = index + 1
        end do
     end do
  end do
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)
  gr_pfftFnArgCopyGridToBuf = 0
end function gr_pfftFnArgCopyGridToBuf
