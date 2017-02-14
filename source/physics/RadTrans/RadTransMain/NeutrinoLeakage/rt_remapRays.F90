!!****if* source/physics/RadTrans/RadTransMain/NeutrinoLeakage/rt_remapRays
!!
!!  NAME 
!!
!!  rt_remapRays
!!
!!  SYNOPSIS
!!
!!  call rt_sampleRays( )
!!
!!  DESCRIPTION 
!!  Adjust centering of leakage rays to be centered on PNS.
!!  Center of the PNS is determined by the maximum density
!!  location.  We subtract off half of the zone spacing so
!!  that the origin is exactly zero if the max dens location
!!  is one of the most central zones.
!!
!!  ARGUMENTS
!!
!!
!!***

subroutine rt_remapRays(nblk,blklst)

#include "Flash.h"
#include "constants.h"

  use rt_data, ONLY :  rt_leakX, rt_leakY, rt_leakZ, rt_leakNumRad, rt_leakNumTht, &
       rt_leakRadii, rt_leakTheta, rt_meshComm, rt_pnsCoord, rt_pnsDens, &
       rt_leakNumPhi, rt_leakPhi
  use RadTrans_data, ONLY : rt_meshMe

  use Grid_interface, ONLY : Grid_releaseBlkPtr, Grid_getBlkPtr, &
       Grid_getBlkBoundBox, Grid_getBlkIndexLimits, Grid_getCellCoords

  implicit none
  include "Flash_mpi.h"

  integer, intent(in) :: nblk
  integer, intent(in) :: blklst(nblk)

  real, pointer :: solnData(:,:,:,:)
  real, dimension(LOW:HIGH,MDIM) :: bndBox
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  integer,dimension(MDIM)  :: dimSize
  real,allocatable, dimension(:) :: xCenter, yCenter, zCenter
  real :: delX, delY, delZ

  integer :: i,j,k,n,blockID
  logical :: gcell = .true.
  real :: localMax, localCoord(MDIM)
  integer :: ierr

  localCoord = 0.
  localMax = 0.

  do n=1,nblk
     blockID = blklst(n)

     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     call Grid_getBlkPtr(blockID,solnData)
     dimSize(:)=blkLimitsGC(HIGH,:)-blkLimitsGC(LOW,:)+1
     if (NDIM > 2)then
        allocate(zCenter(dimSize(KAXIS)))
        call Grid_getCellCoords(KAXIS,blockID,&
             CENTER,gcell,zCenter,dimSize(KAXIS))
        delZ = zCenter(2) - zCenter(1)
     end if
     if (NDIM > 1)then
        allocate(yCenter(dimSize(JAXIS)))
        call Grid_getCellCoords(JAXIS,blockID,&
             CENTER,gcell,yCenter,dimSize(JAXIS))
        delY = yCenter(2) - yCenter(1)
     end if
     allocate(xCenter(dimSize(IAXIS)))
     call Grid_getCellCoords(IAXIS,blockID,&
          CENTER,gcell,xCenter,dimSize(IAXIS))
     delX = xCenter(2) - xCenter(1)

     do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
        do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
           do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
              if (solnData(DENS_VAR,i,j,k) > localMax) then                 
                 localMax = solnData(DENS_VAR,i,j,k)
#if NDIM > 1
                 localCoord(JAXIS) = yCenter(j) - sign(1.0,yCenter(j))*0.5*delY
#if NDIM == 3
                 localCoord(IAXIS) = xCenter(i) - sign(1.0,xCenter(i))*0.5*delX
                 localCoord(KAXIS) = zCenter(k) - sign(1.0,zCenter(k))*0.5*delZ
#endif
#endif
              end if
           end do
        end do
     end do
     
     call Grid_releaseBlkPtr(blockID,solnData)
     deallocate(xCenter)
#if NDIM >1
     deallocate(yCenter)
#if NDIM >2
     deallocate(zCenter)
#endif
#endif
  end do !blkLst

  call MPI_Allreduce(localMax,rt_pnsDens, 1, FLASH_REAL, MPI_MAX, &
                     rt_meshComm, ierr)
  if (rt_pnsDens /= localMax) &
       localCoord = 0.
  call MPI_Allreduce(localCoord, rt_pnsCoord, MDIM, FLASH_REAL, MPI_SUM, &
                     rt_meshComm, ierr)

  do k=1,rt_leakNumPhi
     do j=1,rt_leakNumTht
        do i=1,rt_leakNumRad
#if NDIM > 1
           rt_leakY(i,j,k) = rt_pnsCoord(JAXIS) + rt_leakRadii(i)*cos(rt_leakTheta(j))
#if NDIM == 3
           rt_leakX(i,j,k) = rt_pnsCoord(IAXIS) + rt_leakRadii(i)*sin(rt_leakTheta(j))*cos(rt_leakPhi(k))
           rt_leakZ(i,j,k) = rt_pnsCoord(KAXIS) + rt_leakRadii(i)*sin(rt_leakTheta(j))*sin(rt_leakPhi(k))
#endif
#endif
        end do
     end do
  end do

  return
end subroutine rt_remapRays
