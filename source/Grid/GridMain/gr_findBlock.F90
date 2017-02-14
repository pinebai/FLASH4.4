!!****if* source/Grid/GridMain/gr_findBlock
!!
!! NAME
!!
!!  gr_ptFindBlock
!!
!! SYNOPSIS
!!
!!  call gr_findBlock(integer(in)    :: blkList(blkCount),
!!                    integer(in)    :: blkCount,
!!                    real   (in)    :: pos(MDIM),
!!                    integer(INOUT) :: blockID)
!!
!! DESCRIPTION
!!
!!   Given a point in the domain, this routine finds if the
!!   point lies within one of the blocks on this process. If
!!   such a block is found, its ID is returned, otherwise
!!   blockID is set to NONEXISTENT
!!
!! ARGUMENTS
!!
!!   blkList  : The list of blocks to be examined
!!   blkCount : the number of blocks in the list
!!   pos      : the coordinates of the point
!!   blockID  : the identity of the block if found, NONEXISTENT otherwise
!!              NONEXISTENT is a constant defined in constants.h
!!
!!***
subroutine gr_findBlock(blkList,blkCount,pos,blockID)

#include "constants.h"
#include "Flash.h"

  use Grid_interface, ONLY : Grid_getBlkBoundBox, Grid_outsideBoundBox
  implicit none

  integer,intent(IN) :: blkCount
  integer,dimension(blkCount),intent(IN) :: blkList
  real,dimension(MDIM),intent(IN) :: pos
  integer,intent(INOUT) :: blockID


  real,dimension(LOW:HIGH,MDIM) :: bndBox

  logical :: found,outSide
  integer,dimension(MDIM) :: Negh

  integer :: i,j,k, proc

  found = .false.
  j=0
  do while((.not.found).and.(j<blkCount))
     j=j+1
     blockID=blkList(j)
     call Grid_getBlkBoundBox(blockID,bndBox)
     call Grid_outsideBoundBox(pos,bndBox,outside,Negh)
     found=.not.outSide
  end do
  if(.not.found)then
     blockID=NONEXISTENT
  end if

  return
end subroutine gr_findBlock
