!!****if* source/physics/TreeRay/TreeRayMain/TreeRay_bhSelfContrib
!!
!! NAME
!!
!!  TreeRay_bhSelfContrib
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

subroutine TreeRay_bhSelfContrib(node, cellsize, blockno, point, blkLimits, solnData)
  use TreeRay_data, ONLY : tr_bhIM, tr_bhVolRays, tr_bhMassRays, tr_bhUseTreeRay
  !use tr_osInterface, ONLY : tr_osSelfContrib
  use tr_odInterface, ONLY : tr_odSelfContrib
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
  real, dimension(:), intent(IN) :: node
  real, dimension(MDIM), intent(IN) :: cellsize
  integer, intent(IN) :: blockno
  integer, dimension(MDIM), intent(IN) :: point
  integer, dimension(2,MDIM)   :: blkLimits
  real, DIMENSION(:,:,:,:), POINTER, intent(OUT) :: solnData
  integer :: ii, jj, kk

  if (.not. tr_bhUseTreeRay) return

  ! coordinates in the *_rays arrays
  ii = point(IAXIS) - blkLimits(LOW,IAXIS) + 1
  jj = point(JAXIS) - blkLimits(LOW,JAXIS) + 1
  kk = point(KAXIS) - blkLimits(LOW,KAXIS) + 1

  !call tr_osSelfContrib(node, ii,jj, kk)
  call tr_odSelfContrib(node, ii,jj, kk)


  return
end subroutine TreeRay_bhSelfContrib


