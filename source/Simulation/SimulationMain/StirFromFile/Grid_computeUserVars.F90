!!****if* source/Simulation/SimulationMain/StirFromFile/Grid_computeUserVars
!!
!! NAME
!!  Grid_computeUserVars
!!
!! SYNOPSIS
!!  Grid_computeUserVars()
!!
!! DESCRIPTION
!!
!!  The divergence of the velocity (and random force) and
!!  the magnitude of the vorticity (and random force) 
!!  are computed here.
!!
!!  by Christoph Federrath, 2013
!!
!!***

subroutine Grid_computeUserVars()

  use Grid_interface, ONLY : Grid_getLocalNumBlks, Grid_getBlkIndexLimits, &
                             Grid_getCellCoords, Grid_getBlkPtr, Grid_releaseBlkPtr, &
                             Grid_fillGuardCells

  implicit none

#include "Flash.h"
#include "constants.h"

  integer :: blockID, i, j, k
  integer :: xsize, ysize, zsize, localNumBlocks
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC

  real, DIMENSION(:,:,:,:), POINTER :: solnData

#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_IHI_GC) :: xCenter
  real, dimension(GRID_JHI_GC) :: yCenter
  real, dimension(GRID_KHI_GC) :: zCenter
#else
  real, allocatable, dimension(:) :: xCenter, yCenter, zCenter
  integer :: istat
#endif

  logical, parameter :: Debug = .false.

  if (Debug) print *, 'Grid_computeUserVars:  entering ...'

  ! need to update guard cell info, because we use data in the guard cells below
  call Grid_fillGuardCells(CENTER, ALLDIR)

  call Grid_getLocalNumBlks(localNumBlocks)

  do blockID = 1, localNumBlocks

     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     xsize = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
     ysize = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
     zsize = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

#ifndef FIXEDBLOCKSIZE
     allocate(xCenter(xsize),stat=istat)
     allocate(yCenter(ysize),stat=istat)
     allocate(zCenter(zsize),stat=istat)
#endif

     call Grid_getCellCoords(IAXIS, blockID, CENTER, .true., xCenter, xsize)
     call Grid_getCellCoords(JAXIS, blockID, CENTER, .true., yCenter, ysize)
     call Grid_getCellCoords(KAXIS, blockID, CENTER, .true., zCenter, zsize)

     if (Debug) print *, 'Grid_computeUserVars:  Now calculating vortocity, divergence, etc ...'

     ! get a pointer to the current block of data
     call Grid_getBlkPtr(blockID, solnData)

     do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)

#ifdef DVVL_VAR
              solnData(DVVL_VAR,i,j,k) = &
                (solnData(VELX_VAR,i+1,j,  k  )-solnData(VELX_VAR,i-1,j  ,k  ))/(xCenter(i+1)-xCenter(i-1)) + &
                (solnData(VELY_VAR,i  ,j+1,k  )-solnData(VELY_VAR,i  ,j-1,k  ))/(yCenter(j+1)-yCenter(j-1))
#if NDIM == 3
              solnData(DVVL_VAR,i,j,k) = solnData(DVVL_VAR,i,j,k) + &
                (solnData(VELZ_VAR,i  ,j  ,k+1)-solnData(VELZ_VAR,i  ,j  ,k-1))/(zCenter(k+1)-zCenter(k-1))
#endif
#endif
#ifdef DVRF_VAR
              solnData(DVRF_VAR,i,j,k) = &
                (solnData(ACCX_VAR,i+1,j,  k  )-solnData(ACCX_VAR,i-1,j  ,k  ))/(xCenter(i+1)-xCenter(i-1)) + &
                (solnData(ACCY_VAR,i  ,j+1,k  )-solnData(ACCY_VAR,i  ,j-1,k  ))/(yCenter(j+1)-yCenter(j-1))
#if NDIM == 3
              solnData(DVRF_VAR,i,j,k) = solnData(DVRF_VAR,i,j,k) + &
                (solnData(ACCZ_VAR,i  ,j  ,k+1)-solnData(ACCZ_VAR,i  ,j  ,k-1))/(zCenter(k+1)-zCenter(k-1))
#endif
#endif
#ifdef MVRT_VAR
              solnData(MVRT_VAR,i,j,k) = &
                ( (solnData(VELY_VAR,i+1,j,  k  )-solnData(VELY_VAR,i-1,j  ,k  ))/(xCenter(i+1)-xCenter(i-1)) - &
                  (solnData(VELX_VAR,i  ,j+1,k  )-solnData(VELX_VAR,i  ,j-1,k  ))/(yCenter(j+1)-yCenter(j-1)) )**2
#if NDIM == 3
              solnData(MVRT_VAR,i,j,k) = solnData(MVRT_VAR,i,j,k) + &
                ( (solnData(VELX_VAR,i  ,j,  k+1)-solnData(VELX_VAR,i  ,j  ,k-1))/(zCenter(k+1)-zCenter(k-1)) - &
                  (solnData(VELZ_VAR,i+1,j,  k  )-solnData(VELZ_VAR,i-1,j  ,k  ))/(xCenter(i+1)-xCenter(i-1)) )**2

              solnData(MVRT_VAR,i,j,k) = solnData(MVRT_VAR,i,j,k) + &
                ( (solnData(VELZ_VAR,i  ,j+1,k  )-solnData(VELZ_VAR,i  ,j-1,k  ))/(yCenter(j+1)-yCenter(j-1)) - &
                  (solnData(VELY_VAR,i  ,j,  k+1)-solnData(VELY_VAR,i  ,j  ,k-1))/(zCenter(k+1)-zCenter(k-1)) )**2
#endif
              solnData(MVRT_VAR,i,j,k) = sqrt(solnData(MVRT_VAR,i,j,k))
#endif
#ifdef RTRF_VAR
              solnData(RTRF_VAR,i,j,k) = &
                ( (solnData(ACCY_VAR,i+1,j,  k  )-solnData(ACCY_VAR,i-1,j  ,k  ))/(xCenter(i+1)-xCenter(i-1)) - &
                  (solnData(ACCX_VAR,i  ,j+1,k  )-solnData(ACCX_VAR,i  ,j-1,k  ))/(yCenter(j+1)-yCenter(j-1)) )**2
#if NDIM == 3
              solnData(RTRF_VAR,i,j,k) = solnData(RTRF_VAR,i,j,k) + &
                ( (solnData(ACCX_VAR,i  ,j,  k+1)-solnData(ACCX_VAR,i  ,j  ,k-1))/(zCenter(k+1)-zCenter(k-1)) - &
                  (solnData(ACCZ_VAR,i+1,j,  k  )-solnData(ACCZ_VAR,i-1,j  ,k  ))/(xCenter(i+1)-xCenter(i-1)) )**2

              solnData(RTRF_VAR,i,j,k) = solnData(RTRF_VAR,i,j,k) + &
                ( (solnData(ACCZ_VAR,i  ,j+1,k  )-solnData(ACCZ_VAR,i  ,j-1,k  ))/(yCenter(j+1)-yCenter(j-1)) - &
                  (solnData(ACCY_VAR,i  ,j,  k+1)-solnData(ACCY_VAR,i  ,j  ,k-1))/(zCenter(k+1)-zCenter(k-1)) )**2
#endif
              solnData(RTRF_VAR,i,j,k) = sqrt(solnData(RTRF_VAR,i,j,k))
#endif
           enddo
        enddo
     enddo

     call Grid_releaseBlkPtr(blockID, solnData)

#ifndef FIXEDBLOCKSIZE
     deallocate(xCenter)
     deallocate(yCenter)
     deallocate(zCenter)
#endif

  end do ! loop over blocks

  if (Debug) print *, 'Grid_computeUserVars:  exiting ...'

end subroutine Grid_computeUserVars
