!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhParentContrib
!!
!! NAME
!!
!!  gr_bhParentContrib
!!
!!
!! SYNOPSIS
!!
!!   gr_bhParentContrib(
!!          integer,intent(in) :: block,
!!          integer,intent(in) :: tr,
!!          integer,intent(in) :: cpu,
!!          real, intent(INOUT), dimension(:,:,:,:), pointer :: solnData
!!          )
!!
!! DESCRIPTION
!!
!!   Computes the contribution of a parent block at a specific point.
!!
!! ARGUMENTS
!!
!!  block   : ID of a block where the contribution is calculated
!!  tr      : ID of a block/tree which contributes
!!  cpu     : cpu of a block/tree which contributes
!!  Phi     : 3D array containing potention, contribution added to it
!!
!!***

subroutine gr_bhParentContrib(block, tr, cpu, solnData, blkLimits)

  use gr_bhLocalInterface, ONLY : gr_bhNodeContrib, gr_bhPeriodicDr
  use Driver_interface, ONLY : Driver_abortFlash
  use gr_bhData, ONLY : gr_bhLocCoords, gr_bhTreeNodeSize, &
    gr_bhTreeParentTree, gr_bhTreeLrefine, GR_TREE_IX, GR_TREE_IY, &
    GR_TREE_IZ, GR_TREE_NSIZE

  implicit none
#include "constants.h"
#include "Flash.h"
#include "FortranLangFeatures.fh"

  integer,intent(in) :: block, tr, cpu
  integer :: ii, jj, kk, istat
  integer, dimension(MDIM) :: point
  real, dimension(MDIM+2) :: dr
  real, dimension(:,:,:,:), POINTER_INTENT_IN :: solnData
  integer, dimension(2,MDIM), intent(IN)   :: blkLimits
  real  ::  dist2
  real, allocatable  :: node(:)

  ! block is far enough to be added
  allocate(node(GR_TREE_NSIZE), stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate node in gr_bhParentContrib.F90")
  do kk = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS) ! this block
    point(KAXIS) = kk
    do jj = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
      point(JAXIS) = jj
      do ii = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
        point(IAXIS) = ii
 
        ! compute square of distance
        dr(IAXIS) = gr_bhTreeParentTree(GR_TREE_IX, tr, cpu) &
        &         - gr_bhLocCoords(point(IAXIS)-blkLimits(LOW,IAXIS)+1, IAXIS, block)
        dr(JAXIS) = gr_bhTreeParentTree(GR_TREE_IY, tr, cpu) &
        &         - gr_bhLocCoords(point(JAXIS)-blkLimits(LOW,JAXIS)+1, JAXIS, block)
        dr(KAXIS) = gr_bhTreeParentTree(GR_TREE_IZ, tr, cpu) &
        &         - gr_bhLocCoords(point(KAXIS)-blkLimits(LOW,KAXIS)+1, KAXIS, block)

        ! correct for periodoc conditions
        call gr_bhPeriodicDr(dr)

        dist2  = dr(IAXIS)*dr(IAXIS) + dr(JAXIS)*dr(JAXIS) + dr(KAXIS)*dr(KAXIS)
        dr(MDIM+1) = dist2
        dr(MDIM+2) = sqrt(1.0 / (dist2 + 1D-99))
        node(1:GR_TREE_NSIZE) = gr_bhTreeParentTree(1:GR_TREE_NSIZE, tr, cpu)
        call gr_bhNodeContrib(node, 0, gr_bhTreeLrefine(tr,cpu) &
        & , dr, block, point, blkLimits, solnData)

      enddo
    enddo
  enddo
  deallocate(node)

end subroutine gr_bhParentContrib
