!!****if* source/physics/TreeRay/TreeRayMain/TreeRay_bhFinalizeBlock
!!
!! NAME
!!
!!  TreeRay_bhFinalizeBlock
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

subroutine TreeRay_bhFinalizeBlock(blockno, blkLimits, solnData)
  use Grid_interface, ONLY : Grid_getDeltas
  use TreeRay_data, ONLY : tr_nPix, tr_bhNR, tr_nEb, &
    tr_bhMassRays, tr_bhVolRays, tr_bhSrcfRays, &
    tr_bhEradRays, tr_bhRayR, tr_debprint, tr_meshMe, &
    tr_mapEbSoln, tr_nCd, tr_bhCdMaps, tr_BhUseTreeRay

  use tr_bhLocalInterface, ONLY : tr_bhRadToGas, tr_bhFinalizeCell
  use RuntimeParameters_interface, ONLY: RuntimeParameters_get
  use Eos_interface, ONLY : Eos_wrapped
  use gr_bhData, ONLY : gr_bhLocCoords ! debugging only
  use Timers_interface, ONLY : Timers_start, Timers_stop
  implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
#include "Eos.h"

  integer, intent(IN) :: blockno
  integer, dimension(2,MDIM)   :: blkLimits
  real, DIMENSION(:,:,:,:), POINTER :: solnData
  real, dimension(:), pointer :: solnPoint
  integer :: i,j,k, ipix, ii,jj,kk, ir, ieb
  real :: rho_poc, temp_poc, mass_poc, vol_poc, area_poc
  real :: eflux(tr_nEb, 0:tr_nPix-1)
  real :: cdMaps(tr_nCd, 0:tr_nPix-1)
  real :: rho_ray(0:tr_bhNR), erad_ray(tr_nEb, 0:tr_bhNR)
  real :: del(MDIM)
  real :: dt
  logical :: raycheck

  ! DEBUG
  integer :: ins, ith, iph, iii
  real :: ns, theta, phi, zor, weight, mass_tot, vol_tot

  if (.not. tr_bhUseTreeRay) return

  call Timers_start("tr_finalize_block")

  call Driver_getDt(dt)
  call Grid_getDeltas(blockno,del)
  vol_poc = del(IAXIS)*del(JAXIS)*del(KAXIS)
  area_poc = del(IAXIS)*del(IAXIS)

  do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
    do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
      do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)

        ii = i - blkLimits(LOW,IAXIS) + 1
        jj = j - blkLimits(LOW,JAXIS) + 1
        kk = k - blkLimits(LOW,KAXIS) + 1


#ifdef TR_OPTICALDEPTH
        cdMaps = tr_bhCdMaps(:,:,ii,jj,kk)
#endif
        solnPoint => solnData(:,i,j,k)
        call tr_bhFinalizeCell(solnPoint, del, eflux, cdMaps)

      enddo
    enddo
  enddo

  call Timers_stop("tr_finalize_block")

end subroutine TreeRay_bhFinalizeBlock
