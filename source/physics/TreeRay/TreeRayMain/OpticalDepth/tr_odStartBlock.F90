!!****if* source/physics/TreeRay/TreeRayMain/OpticalDepth/tr_odStartBlock
!!
!! NAME
!!
!!  tr_odStartBlock
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

subroutine tr_odStartBlock(blockno, blkLimits, solnData)
  use TreeRay_data, ONLY : tr_bhCdMaps
  use tr_odData, ONLY : tr_odICDTO, tr_odICDH2, tr_odICDCO
  implicit none
#include "constants.h"
#include "Flash.h"
  integer, intent(IN) :: blockno
  integer, dimension(2,MDIM), intent(IN)   :: blkLimits
  real, DIMENSION(:,:,:,:), POINTER :: solnData

  tr_bhCDMaps(tr_odICDTO,:,:,:,:) = 0.0
#if CHEMISTRYNETWORK == 4 || CHEMISTRYNETWORK == 5
  tr_bhCDMaps(tr_odICDH2,:,:,:,:) = 0.0
#endif
#if CHEMISTRYNETWORK == 5
  tr_bhCDMaps(tr_odICDCO,:,:,:,:) = 0.0
#endif

  return
end subroutine tr_odStartBlock

