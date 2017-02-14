!!****if* source/physics/TreeRay/TreeRayMain/TreeRay_bhNodeContrib
!!
!! NAME
!!
!!  TreeRay_bhNodeContrib
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

subroutine TreeRay_bhNodeContrib(node, trLevel, refLevel, dr, blockno, point, blkLimits, solnData)
  use TreeRay_data, ONLY : tr_bhUseTreeRay, tr_ilNTheta, tr_ilNPhi, tr_ilNNS, tr_ilNI, &
    & tr_ilNNS, tr_ilNSSampFacI, tr_intersectList, tr_bhNR, tr_bhMaxDist, tr_bhMassRays, &
    tr_bhVolRays, tr_bhIM, tr_bhIBM, tr_bhIX, tr_bhIY, tr_bhIZ, &
    tr_nFineR, tr_bhRadNodeMapInd, tr_bhRadNodeMapVal, tr_bhTreeLevels, &
    tr_bhMinCellSize, tr_bhNTBLevels, tr_bhRayRadRes, tr_bhTreeNodeSize

  !use tr_osInterface, ONLY : tr_osNodeContrib
  use tr_odInterface, ONLY : tr_odNodeContrib
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
  real, dimension(:), intent(IN) :: node
  integer, intent(IN) :: trLevel, refLevel, blockno
  integer, dimension(MDIM), intent(IN) :: point
  real, dimension(MDIM+2), intent(IN) :: dr
  integer, dimension(2,MDIM), intent(IN)   :: blkLimits
  real, DIMENSION(:,:,:,:), POINTER :: solnData

  integer :: i, j, ii, jj, kk, ir, irf, ins, ith, iph, ipix
  integer :: ndLevel, iNodeSize
  real :: ndSize, ndMass, ndVol, weight, drmag, ns
  real :: theta, phi, zor, tot_weight, tot_w2
  real :: xtest, ytest, ztest, dotprod

  if (.not. tr_bhUseTreeRay) return

  ! check if the distance does not exceed tr_maxdist
  drmag = sqrt(dr(MDIM+1))
  if (drmag > tr_bhMaxDist) return

  ! node size
  ndLevel = trLevel + refLevel
  ndSize = tr_bhTreeNodeSize(ndLevel)

  ! coordinates in the *_rays arrays
  ii = point(IAXIS) - blkLimits(LOW,IAXIS) + 1
  jj = point(JAXIS) - blkLimits(LOW,JAXIS) + 1
  kk = point(KAXIS) - blkLimits(LOW,KAXIS) + 1

  ! node angular size - assumes that the node is cubic
  ns = ndSize * dr(MDIM+2)
  ins = floor(ns*tr_ilNNS*tr_ilNSSampFacI + 0.5)
  if ((ins < 1) .or. (ins > tr_ilNNS)) then
    print *, "Wrong ins: ", ins, tr_ilNNS, ns, floor(ns*tr_ilNNS + 0.5)
    call Driver_abortFlash("TreeRay_bhNodeContrib: incorrect node size index")
  endif

  ! node angular position
  zor = max(-1.0, min(1.0, dr(KAXIS)*dr(MDIM+2)))
  theta = acos(zor)
  phi = mod(atan2(dr(JAXIS), dr(IAXIS))+2*PI, 2*PI)
  ith = floor(theta*tr_ilNTheta/PI+0.5)
  iph = floor(0.5*phi*tr_ilNPhi/PI+0.5)
  if ((ith < 0) .or. (ith > tr_ilNTheta)) then
    print *, "NC Wrong ith: ", ith, tr_ilNTheta, theta, floor(theta*tr_ilNTheta/PI+0.5), dr
    call Driver_abortFlash("TreeRay_bhBotNodeContrib: incorrect theta index")
  endif
  if ((iph < 0) .or. (iph > tr_ilNPhi)) then
    print *, "Wrong iph: ", iph, tr_ilNPhi, phi, floor(0.5*phi*tr_ilNPhi/PI+0.5)
    call Driver_abortFlash("TreeRay_bhBotNodeContrib: incorrect phi index")
  endif

  call tr_odNodeContrib(node, trLevel, refLevel, ii, jj, kk, ins, iph, ith, dr)
  !call tr_osNodeContrib(node, trLevel, refLevel, ii, jj, kk, ins, iph, ith, irf)

  return


end subroutine TreeRay_bhNodeContrib


