!!****if* source/physics/TreeRay/TreeRayMain/OpticalDepth/tr_odAccNode
!!
!! NAME
!!
!!  tr_odAccNode
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

subroutine tr_odAccNode(subnode, accnode)
  use tr_odData, ONLY : tr_odIH2, tr_odICO
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
  real, dimension(:), intent(IN)  :: subnode
  real, dimension(:), intent(OUT) :: accnode

# if CHEMISTRYNETWORK == 4 || CHEMISTRYNETWORK == 5
  accnode(tr_odIH2) = accnode(tr_odIH2) + subnode(tr_odIH2)
#endif
#if CHEMISTRYNETWORK == 5
  accnode(tr_odICO) = accnode(tr_odICO) + subnode(tr_odICO)
#endif

end subroutine tr_odAccNode
