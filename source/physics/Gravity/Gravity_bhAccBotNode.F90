!!****f* source/physics/Gravity/Gravity_bhAccBotNode
!!
!! NAME
!!
!!  Gravity_bhAccBotNode
!!
!!
!! SYNOPSIS
!!
!!   call Gravity_bhAccBotNode(
!!                          integer(in)    :: blockno,
!!                          integer(in)    :: point(MDIM),
!!                          integer(in)    :: blkLimits(2,MDIM),
!!                          real(in)       :: locCoords(:,:,:),
!!                          real,pointer   :: solnData(:,:,:,:),
!!                          real(in)       :: botnode(:),
!!                          real(inout)    :: accnode(:)
!!        )
!!
!! DESCRIPTION
!!
!!  Called during tree build. Adds values of botnode into accnode.
!!
!! ARGUMENTS
!!
!!  blockno     : number of block into which the node belongs
!!  point       : indeces of the cell (botnode) in the block
!!  blkLimits   : limits of indeces in the block
!!  locCoords   : array with coordinates of allblocks on a given cpu
!!  solnData    : solution data from the grid
!!  botnode     : array of the bottom-most (leaf) node of the tree, i.e. a grid cell
!!  accnode     : array of the node into which botnode contribution is added
!!
!!
!!***

subroutine Gravity_bhAccBotNode(blockno, point, blkLimits, locCoords, solnData, botnode, accnode)
  implicit none
#include "constants.h"
#include "Flash.h"
  integer, intent(IN) :: blockno
  integer, dimension(MDIM), intent(IN) :: point
  integer, dimension(2,MDIM), intent(IN)   :: blkLimits
  real, dimension(:,:,:), intent(IN) :: locCoords
  real, DIMENSION(:,:,:,:), POINTER :: solnData
  real, dimension(:), intent(IN) :: botnode
  real, dimension(:), intent(INOUT) :: accnode
  return
end subroutine Gravity_bhAccBotNode
