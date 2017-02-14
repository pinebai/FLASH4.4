!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhMAC
!!
!! NAME
!!
!!  gr_bhMAC
!!
!!
!! SYNOPSIS
!!
!!  logical res = gr_bhMAC(
!!                          logical(in)    :: physMAC,
!!                          real(in)       :: node(:),
!!                          real(in)       :: ndSize2,
!!                          real(in)       :: dr(MDIM+2),
!!                          integer(in)    :: blockno,
!!                          integer(in)    :: point(MDIM),
!!                          integer(in)    :: blkLimits(2,MDIM),
!!                          real,pointer   :: solnData(:,:,:,:)
!!        )
!!
!! DESCRIPTION
!!
!!  Multipole Acceptance Criterion. Determines whether the contribution of the
!!  node to the potential at the point of calculation will have a nacessary
!!  accuracy.
!!
!! ARGUMENTS
!!
!!  physMAC     : if TRUE, MACs of physical modules will be used
!!                if FALSE, only purely geometric (Barnes-Hut) MAC is used
!!  node        : array of the node tested
!!  ndSize2     : square of the physical size of the node (the largest 
!!                extent of the node)
!!  dr          : (1:MDIM) - position vector from the point-of-calculation to the node
!!                (MDIM+1) - square of the magnitude of the position vector
!!                (MDIM+2) - inverted magnitude of the position vector
!!  blockno     : number of block into which the point-of-calculation belongs
!!  point       : indeces of the point-of-calculation in the block
!!  blkLimits   : limits of indeces in the block
!!  solnData    : solution data from the grid
!!
!! RESULT
!!  Returns TRUE if the node is accepted for the calculation of the potential.
!!  Otherwise returns FALSE.
!!
!!
!!***

logical function gr_bhMAC(physMAC, node, ndSize2, drGC, dr, blockno, point, blkLimits, solnData)
  use gr_bhData, ONLY : gr_bhTreeLimangle2i, gr_bhTreeSafeBox, &
    & gr_bhTreeSafeBoxHalf2, gr_bhTreeMyPE, gr_bhUseGravity, gr_bhUseTreeRay, &
    & gr_bhPhysMACComm_step
  use Gravity_interface, ONLY : Gravity_bhMAC
  use TreeRay_interface, ONLY : TreeRay_bhMAC
  implicit none
#include "FortranLangFeatures.fh"
#include "constants.h"
  logical, intent(IN) :: physMAC
  real, dimension(:), intent(IN) :: node
  real, intent(IN) :: ndSize2
  integer, dimension(MDIM), intent(IN) :: point
  real, dimension(MDIM+2), intent(IN) :: dr
  real, dimension(MDIM), intent(IN) :: drGC
  integer, intent(IN) :: blockno
  integer, dimension(2,MDIM), intent(IN)   :: blkLimits
  real, dimension(MDIM) :: drGC2
  real, DIMENSION(:,:,:,:), POINTER_INTENT_IN :: solnData
  logical :: res
  real :: ndSize2I

  res = .true.

  ! check if the point is not inside the safe box
  if (gr_bhTreeSafeBox > 0) then
    ndSize2I = 1.0/ndSize2
    drGC2 = drGC*drGC
    if (    (drGC2(IAXIS)*ndSize2I < gr_bhTreeSafeBoxHalf2) &
    & .and. (drGC2(JAXIS)*ndSize2I < gr_bhTreeSafeBoxHalf2) &
    & .and. (drGC2(KAXIS)*ndSize2I < gr_bhTreeSafeBoxHalf2)) then
      !print *, "MAC: ", drGC2(:), ndSize2I, gr_bhTreeSafeBoxHalf2, drGC2(:)*ndSize2I
      res = .false.
    endif
  endif

  if (.not. physMAC) then
    ! check the opening angle criterion
    if (dr(MDIM+1) < ndSize2*gr_bhTreeLimangle2i) then
      res = .false.
    endif
  else
    if (.not. gr_bhPhysMACComm_step .and. dr(MDIM+1) > ndSize2*gr_bhTreeLimangle2i) then
      ! if geometric criterion was used for communication, accept the node that
      ! fullfills it; otherwise, node may be not present on a given CPU
      res = .true.
    else
      ! MAC of physical modules
      ! res is tested because fortran does not ensure short-circuit .AND. evaluation
      if (gr_bhUseGravity .and. res) &
      & res = res .and. Gravity_bhMAC(node, ndSize2, dr, blockno, point, blkLimits, solnData)
      if (gr_bhUseTreeRay .and. res) &
      & res = res .and. TreeRay_bhMAC(node, ndSize2, dr, blockno, point, blkLimits, solnData)
    endif
  endif ! physMAC

  !if ((gr_bhTreeMyPE == 0) .and. (blockno == 44) &
  !& .and. (point(IAXIS) == 5) .and. (point(JAXIS) == 5).and. (point(KAXIS) == 5)) then
  !  print *, "MAC3: ", res
  !endif

  gr_bhMAC = res

end function gr_bhMAC

