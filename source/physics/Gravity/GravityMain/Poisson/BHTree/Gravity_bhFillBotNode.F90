!!****if* source/physics/Gravity/GravityMain/Poisson/BHTree/Gravity_bhFillBotNode
!!
!! NAME
!!
!!  Gravity_bhFillBotNode
!!
!!
!! SYNOPSIS
!!
!!   call Gravity_bhFillBotNode(
!!                          integer(in)    :: blockno,
!!                          integer(in)    :: point(MDIM),
!!                          integer(in)    :: blkLimits(2,MDIM),
!!                          real,pointer   :: solnData(:,:,:,:),
!!                          real(inout)    :: botnode(:)
!!        )
!!
!! DESCRIPTION
!!
!!  Called during tree build. Fills the bottom-most (leaf) node 
!!  (corresponding to a single grid cell) with values from the grid.
!!  Recently empty, because evrything is done in gr_bhFillBotNode.
!!
!! ARGUMENTS
!!
!!  blockno     : number of block into which the node belongs
!!  point       : indeces of the cell (botnode) in the block
!!  blkLimits   : limits of indeces in the block
!!  solnData    : solution data from the grid
!!  botnode     : array of the bottom-most (leaf) node of the tree
!!                which is filled by values from the grid
!!
!!***

subroutine Gravity_bhFillBotNode(blockno, point, blkLimits, solnData, botnode)
  use Gravity_data, ONLY : useGravity !, grv_bhIB2, grv_bhNODE5
  use Grid_interface, ONLY : Grid_getDeltas
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
  integer, intent(IN) :: blockno
  integer, dimension(MDIM), intent(IN) :: point
  integer, dimension(2,MDIM)   :: blkLimits
  real, DIMENSION(:,:,:,:), POINTER :: solnData
  real, dimension(:), intent(INOUT) :: botnode

  if (.not. useGravity) return


  ! NOT USED: in recent version, botnodes do not have a field for 
  ! storing the second order B2 moment / minimum distance
  ! B2 of bottom nodes is calculated in AccBotNode and added directly 
  ! to B2 of the second bottom most level
  !if (grv_bhNODE5 /= grv_bhN5_NONE) then
  !  call Grid_getDeltas(blockno, del)
  !  botnode(grv_bhIB2) = (1.0/12.0) * botnode(grv_bhIM) &
  !  & * (del(IAXIS)**2 + del(JAXIS)**2 + del(KAXIS)**2)
  !  print *, "FBN: ", botnode
  !endif

  return
end subroutine Gravity_bhFillBotNode
