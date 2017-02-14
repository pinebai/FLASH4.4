!!****if* source/physics/Gravity/GravityMain/Poisson/BHTree/Gravity_bhInitFieldVar
!!
!! NAME
!!
!!  Gravity_bhInitFieldVar
!!
!!
!! SYNOPSIS
!!
!!  call Gravity_bhInitFieldVar(
!!                 integer(in) :: gpotVar
!!  )
!!
!! DESCRIPTION
!!
!!  Called before the tree walk. Initializes fiels variables before the potential
!!  is calculated. Sets index of the grid variable where the gravitational 
!!  potential is stored (passed from Grid_solvePoisson - ipotvar. Calculate
!!  gravitational acceleration at each grid cell (inverted maximum allowed error
!!  in acceleration stored in field variable ACEI) - needed for normalization of 
!!  the error.
!!
!!
!! ARGUMENTS
!!
!!  gpotVar : index of the grid variable where the gravitational potential is
!!            stored. Passed from Grid_solvePoisson.
!!
!!
!!***

subroutine Gravity_bhInitFieldVar(gpotVar)
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr &
  & , Grid_getListOfBlocks, Grid_getBlkPhysicalSize, Grid_getBlkIndexLimits
  !use Grid_interface, ONLY : Grid_fillGuardCells
  use Gravity_data, ONLY : useGravity, grv_bhUseRelAccErr, grv_bhAccErr, &
    grv_defaultGpotVar
  use Driver_interface, ONLY : Driver_getNStep

  implicit none
#include "constants.h"
#include "Flash.h"
  integer, intent(IN) :: gpotVar
  integer :: blockCount, blockID
  integer :: blockList(MAXBLOCKS)
  real, DIMENSION(:,:,:,:), POINTER :: solnData
  real :: blockSize(MDIM)
  integer, dimension(2,MDIM)   :: blkLimits,blkLimitsGC
  real :: delx_inv, dely_inv, delz_inv
  integer :: ii, jj, kk, nstep

  if (.not. useGravity) return

  ! set grv_defaultGpotVar for this step
  grv_defaultGpotVar = gpotVar

  ! derivatives of GPOT are calculated for ACEI array
  ! it may be good idea to fill guard cells
  ! call Grid_fillGuardCells(CENTER, ALLDIR)

  call Grid_getListOfBlocks(LEAF,blockList,blockCount) ! get list of LEAF blocks
  do blockID = 1, blockCount
    call Grid_getBlkIndexLimits(blockList(blockID),blkLimits,blkLimitsGC) 
    call Grid_getBlkPtr(blockList(blockID),solnData,CENTER)
    call Grid_getBlkPhysicalSize(blockList(blockID), blockSize)

    ! calculate the gravitational acceleration
    ! determine the maximum allowed error for the node contribution

    !if (grav_useRelAccErr) then
    ! temporarily for all criteria, because it is used for accuracy tests
      delx_inv = (blkLimits(HIGH,IAXIS) - blkLimits(LOW,IAXIS) + 1) / blockSize(IAXIS)
      dely_inv = (blkLimits(HIGH,JAXIS) - blkLimits(LOW,JAXIS) + 1) / blockSize(JAXIS)
      delz_inv = (blkLimits(HIGH,KAXIS) - blkLimits(LOW,KAXIS) + 1) / blockSize(KAXIS)
      do kk = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS) ! this block
        do jj = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
          do ii = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
            solnData(ACEI_VAR,ii,jj,kk) = 2.0 / (sqrt( &
            & ((solnData(grv_defaultGpotVar,ii-1,jj,kk) - solnData(grv_defaultGpotVar,ii+1,jj,kk))*delx_inv)**2 + &
            & ((solnData(grv_defaultGpotVar,ii,jj-1,kk) - solnData(grv_defaultGpotVar,ii,jj+1,kk))*dely_inv)**2 + &
            & ((solnData(grv_defaultGpotVar,ii,jj,kk-1) - solnData(grv_defaultGpotVar,ii,jj,kk+1))*delz_inv)**2) &
            & * grv_bhAccErr + 1d-99)
          enddo
        enddo
      enddo
    !endif

    ! set gpot/gac[xyz] field arrays to zero
#ifdef GRAV_TREE_ACC
    solnData(GACX_VAR, :, :, :) = 0.0
    solnData(GACY_VAR, :, :, :) = 0.0
    solnData(GACZ_VAR, :, :, :) = 0.0
#endif
    solnData(grv_defaultGpotVar, :, :, :) = 0.0
    call Grid_releaseBlkPtr(blockList(blockID),solnData)
  enddo
  
  return
end subroutine Gravity_bhInitFieldVar



