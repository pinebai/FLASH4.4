!!****if* source/Grid/GridMain/Chombo/gr_findBlock
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

#include "constants.h"
#include "Flash.h"

subroutine gr_findBlock(blkList,blkCount,pos,blockID)
  use Driver_interface, ONLY : Driver_abortFlash
  implicit none
  integer,intent(IN) :: blkCount
  integer,dimension(blkCount),intent(IN) :: blkList
  real,dimension(MDIM),intent(IN) :: pos
  integer,intent(INOUT) :: blockID

  call Driver_abortFlash("Not yet implemented for Chombo")
end subroutine gr_findBlock
