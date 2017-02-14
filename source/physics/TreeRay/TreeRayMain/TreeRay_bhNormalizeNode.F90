!!****if* source/physics/TreeRay/TreeRayMain/TreeRay_bhNormalizeNode
!!
!! NAME
!!
!!  TreeRay_bhNormalizeNode
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

subroutine TreeRay_bhNormalizeNode(smr, node)
  use TreeRay_data, ONLY : tr_BhUseTreeRay
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
  real, dimension(MDIM), intent(IN) :: smr
  real, dimension(:), intent(OUT) :: node

  if (.not. tr_bhUseTreeRay) return

  !if (node(tr_bhISrc) > 0.0) then
  !  print *, "NN: node_srcf: ", node
  !endif
  !if (node(tr_bhISrc) > 0.0) then
  !  node(tr_bhISrcX) = node(tr_bhISrcX) / (node(tr_bhISrc) + 1d-99)
  !  node(tr_bhISrcY) = node(tr_bhISrcY) / (node(tr_bhISrc) + 1d-99)
  !  node(tr_bhISrcZ) = node(tr_bhISrcZ) / (node(tr_bhISrc) + 1d-99)
  !else
  !  node(tr_bhISrcX:tr_bhISrcZ) = node(tr_bhIX:tr_bhIZ)
  !endif

  return
end subroutine TreeRay_bhNormalizeNode

