!!****if* source/physics/TreeRay/TreeRayMain/OpticalDepth/tr_odFillBotNode
!!
!! NAME
!!
!!  tr_odFillBotNode
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

subroutine tr_odFillBotNode(blockno, point, blkLimits, solnData, botnode)
  use tr_odData, ONLY : tr_odIBH2, tr_odIBCO
  use Grid_interface, ONLY : Grid_getSingleCellVol
  implicit none
#include "constants.h"
#include "Flash.h"

  integer, intent(IN) :: blockno
  integer, dimension(MDIM), intent(IN) :: point
  integer, dimension(2,MDIM)   :: blkLimits
  real, DIMENSION(:,:,:,:), POINTER, intent(IN) :: solnData
  real, dimension(:), intent(OUT) :: botnode
  real :: dvol

  call Grid_getSingleCellVol(blockno, INTERIOR, point, dvol)
# if CHEMISTRYNETWORK == 4 || CHEMISTRYNETWORK == 5
  !! see Chemistry.F90 for required order of species
  botnode(tr_odIBH2) = solnData(IH2_SPEC,point(IAXIS),point(JAXIS), point(KAXIS)) &
  & * solnData(DENS_VAR,point(IAXIS),point(JAXIS), point(KAXIS)) * dvol
#endif
#if CHEMISTRYNETWORK == 5
  botnode(tr_odIBCO) = solnData(ICO_SPEC,point(IAXIS),point(JAXIS), point(KAXIS)) &
  & * solnData(DENS_VAR,point(IAXIS),point(JAXIS), point(KAXIS)) * dvol
#endif

  return
end subroutine tr_odFillBotNode

