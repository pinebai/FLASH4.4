!!****if* source/Simulation/SimulationMain/StirTurb/Grid_computeUserVars
!!
!! NAME
!!  Grid_computeUserVars
!!
!!
!! SYNOPSIS
!!
!!  call Grid_computeUserVars() 
!!
!!  
!! DESCRIPTION 
!!  
!!  This routine is always a stub unless the user adds this file to their simulations
!!  directory.  Routine allows users to compute own variables.  Some examples are
!!  vorticity, conductivity etc. 
!!
!!  For vorticity, the components of the curl depend on velocity components
!!  transverse to a given direction (eg, x and y for the z-component of the curl).
!!  In the code, these components are defined in cut planes transverse to a given
!!  component of the curl (eg, XY plane for the z-component of the curl), and
!!  then looped over the transverse direction (eg, z for the z-component).
!!  Therefore, there are two definitions for each velocity component,
!!  one for each of two planes (eg, vz defined in the XZ plane and vz defined
!!  in the YZ plane -- vz_xy and vz_yz in the code).
!!
!!
!! ARGUMENTS 
!!  
!!  none
!!
!!
!!
!! NOTES
!!  
!!***

subroutine Grid_computeUserVars()
  use Grid_interface, ONLY : Grid_getLocalNumBlks, &
    Grid_getBlkIndexLimits, Grid_getCellCoords, Grid_getPlaneData, &
    Grid_putBlkData

  implicit none

#include "Flash.h"
#include "constants.h"
  
  integer :: blockID, i, j, k
  integer :: xsize, ysize, zsize, localNumBlocks
  integer,dimension(2) ::   planeSize
  integer,dimension(MDIM) :: startPos, dataSize
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
#ifdef FIXEDBLOCKSIZE

  real, dimension(GRID_IHI_GC) :: xCenter
  real, dimension(GRID_JHI_GC) :: yCenter
  real, dimension(GRID_KHI_GC) :: zCenter

! Cut planes of components of velocity field used in computing gradient

  real, dimension(GRID_JHI_GC, GRID_KHI_GC) :: vy_yz, vz_yz, xvortplane
  real, dimension(GRID_IHI_GC, GRID_KHI_GC) :: vx_xz, vz_xz, yvortplane
  real, dimension(GRID_IHI_GC, GRID_JHI_GC) :: vx_xy, vy_xy, zvortplane
 
  real, dimension(GRID_IHI_GC, GRID_JHI_GC, GRID_KHI_GC) :: &
       mvort, xvort, yvort, zvort
#else
  real,allocatable,dimension(:)::xCenter,yCenter,zCenter
  real,allocatable,dimension(:,:)::vy_yz, vz_yz, xvortplane,&
                                   vx_xz, vz_xz, yvortplane,&
                                   vx_xy, vy_xy, zvortplane
  real,allocatable,dimension(:,:,:)::mvort, xvort, yvort, zvort
#endif


  call Grid_getLocalNumBlks(localNumBlocks)

  do blockID=1,localNumBlocks

     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     xsize = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
     ysize = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
     zsize = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

#ifndef FIXEDBLOCKSIZE
     allocate(xCenter(xsize))
     allocate(yCenter(ysize))
     allocate(zCenter(zsize))
     allocate(vy_yz(ysize,zsize))
     allocate(vz_yz(ysize,zsize))
     allocate(xvortplane(ysize,zsize))
     allocate(vx_xz(xsize,zsize))
     allocate(vz_xz(xsize,zsize))
     allocate(yvortplane(xsize,zsize))
     allocate(vx_xy(xsize,ysize))
     allocate(vy_xy(xsize,ysize))
     allocate(zvortplane(xsize,ysize))
     allocate(mvort(xsize,ysize,zsize))
     allocate(xvort(xsize,ysize,zsize))
     allocate(yvort(xsize,ysize,zsize))
     allocate(zvort(xsize,ysize,zsize))
#endif

     call Grid_getCellCoords(IAXIS, blockID, CENTER, .true., xCenter, xsize)
     call Grid_getCellCoords(JAXIS, blockID, CENTER, .true., yCenter, ysize)
     call Grid_getCellCoords(KAXIS, blockID, CENTER, .true., zCenter, zsize)
     
     !---------------------------------------------------------------------
     ! x-vorticity -- computed over YZ planes; loop over transverse x direction
     !---------------------------------------------------------------------
     xvortplane = 0.0

#if NDIM == 3
     do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
        startPos(1) = i
        startPos(2) = 1
        startPos(3) = 1
        planeSize(1) = blkLimitsGC(HIGH,JAXIS)
        planeSize(2) = blkLimitsGC(HIGH,KAXIS)
        
        call Grid_getPlaneData(blockID, CENTER, VELY_VAR, EXTERIOR, YZPLANE, &
             startPos, vy_yz, planeSize)
        call Grid_getPlaneData(blockID, CENTER, VELZ_VAR, EXTERIOR, YZPLANE, &
             startPos, vz_yz, planeSize)
        do j=blkLimitsGC(LOW,JAXIS)+1,blkLimitsGC(HIGH,JAXIS)-1
           do k=blkLimitsGC(LOW,KAXIS)+1,blkLimitsGC(HIGH,KAXIS)-1
              xvortplane (j,k) = (vz_yz(j+1,k)-vz_yz(j-1,k))/&
                   (yCenter(j+1)-yCenter(j-1)) &
                   -(vy_yz(j,k+1)-vy_yz(j,k-1))/&
                   (zCenter(k+1)-zCenter(k-1))
              
           enddo
        enddo
        
 
        xvort(i,:,:) = xvortplane (:,:)
        
     enddo
     
#endif     
     
     !--------------------------------------------------------------------
     ! y-vorticity -- computed over XZ planes; loop over transverse y direction
     !------------------------------------------------------------------------
     
     yvortplane = 0.0

#if NDIM ==3
     do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        startPos(1) = 1
        startPos(2) = j
        startPos(3) = 1
        planeSize(1) = blkLimitsGC(HIGH,IAXIS)
        planeSize(2) = blkLimitsGC(HIGH,KAXIS)
        
        call Grid_getPlaneData(blockID, CENTER, VELX_VAR, EXTERIOR, XZPLANE, &
             startPos, vx_xz, planeSize)
        call Grid_getPlaneData(blockID, CENTER, VELZ_VAR, EXTERIOR, XZPLANE, &
             startPos, vz_xz, planeSize)
        
        do k=blkLimitsGC(LOW,KAXIS)+1,blkLimitsGC(HIGH,KAXIS)-1
           do i=blkLimitsGC(LOW,IAXIS)+1,blkLimitsGC(HIGH,IAXIS)-1
              yvortplane (i,k) = (vx_xz(i,k+1)-vx_xz(i,k-1))/(zCenter(k+1)-zCenter(k-1)) &
                   - (vz_xz(i+1,k)-vz_xz(i-1,k))/(xCenter(i+1)-xCenter(i-1)) 
           enddo
        enddo
        
        
        yvort(:,j,:) = yvortplane (:,:)
        
     enddo
#endif     
     
     
     !-----------------------------------------------------------------------
     ! z-vorticity -- computed over XY planes; loop over transverse z direction
     !-----------------------------------------------------------------------
     
     zvortplane  = 0.0
     do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
        startPos(1) = 1
        startPos(2) = 1
        startPos(3) = k
        planeSize(1) = blkLimitsGC(HIGH,IAXIS)
        planeSize(2) = blkLimitsGC(HIGH,JAXIS)
        
        call Grid_getPlaneData(blockID, CENTER, VELX_VAR, EXTERIOR, XYPLANE, &
             startPos, vx_xy, planeSize)
        call Grid_getPlaneData(blockID, CENTER, VELY_VAR, EXTERIOR, XYPLANE, &
             startPos, vy_xy, planeSize)
        
        
        do j=blkLimitsGC(LOW,JAXIS)+1,blkLimitsGC(HIGH,JAXIS)-1
           do i=blkLimitsGC(LOW,IAXIS)+1,blkLimitsGC(HIGH,IAXIS)-1
              zvortplane (i,j) = (vy_xy(i+1,j)-vy_xy(i-1,j))/(xCenter(i+1)-xCenter(i-1)) &
                   -(vx_xy(i,j+1)-vx_xy(i,j-1))/(yCenter(j+1)-yCenter(j-1))
           enddo
        enddo
        
        zvort(:,:,k) = zvortplane (:,:)
        
     enddo
     
     do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)
        do j = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
           do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
              
              mvort (i, j, k) = sqrt ( xvort (i, j, k)**2. + yvort (i, j, k)**2. &
                   +  zvort (i, j, k)**2. )
           end do
        end do
     end do
     
     datasize(1) = blkLimitsGC(HIGH,IAXIS)
     datasize(2) = blkLimitsGC(HIGH,JAXIS)
     datasize(3) = blkLimitsGC(HIGH,KAXIS)
     
     startPos=1
     call Grid_putBlkData(blockID, SCRATCH_CTR, MVRT_SCRATCH_CENTER_VAR, EXTERIOR,&
          startPos,mvort, datasize) 
     
  
  
#ifndef FIXEDBLOCKSIZE
     deallocate(xCenter)
     deallocate(yCenter)
     deallocate(zCenter)
     deallocate(vy_yz)
     deallocate(vz_yz)
     deallocate(xvortplane)
     deallocate(vx_xz)
     deallocate(vz_xz)
     deallocate(yvortplane)
     deallocate(vx_xy)
     deallocate(vy_xy)
     deallocate(zvortplane)
     deallocate(mvort)
     deallocate(xvort)
     deallocate(yvort)
     deallocate(zvort)
#endif


  end do

end subroutine Grid_computeUserVars
