!!****if* source/physics/TreeRay/TreeRayMain/OpticalDepth/tr_odAccBotNode
!!
!! NAME
!!
!!  tr_odAccBotNode
!!
!!
!! SYNOPSIS
!!
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!
!! RESULT
!!
!!
!!***

subroutine tr_odAccBotNode(blockno, point, blkLimits, locCoords, solnData, &
   & botnode, accnode)
  use tr_odData, ONLY : tr_odIH2, tr_odIBH2, tr_odICO, tr_odIBCO
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  integer, intent(IN) :: blockno
  integer, dimension(MDIM), intent(IN) :: point
  integer, dimension(2,MDIM)   :: blkLimits
  real, dimension(:,:,:), intent(IN) :: locCoords
  real, DIMENSION(:,:,:,:), POINTER, intent(IN) :: solnData
  real, dimension(:), intent(IN) :: botnode
  real, dimension(:), intent(INOUT) :: accnode

# if CHEMISTRYNETWORK == 4 || CHEMISTRYNETWORK == 5
  accnode(tr_odIH2) = accnode(tr_odIH2) + botnode(tr_odIBH2)
#endif
#if CHEMISTRYNETWORK == 5
  accnode(tr_odICO) = accnode(tr_odICO) + botnode(tr_odIBCO)
#endif

end subroutine tr_odAccBotNode
