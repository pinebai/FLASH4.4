!!****if* source/Grid/GridSolvers/IsoBndMultipole/gr_isoSumLocal
!!
!! NAME
!!
!!  gr_isoSumLocal
!!
!! 
!! SYNOPSIS
!!
!!  call gr_isoSumLocal(real(INOUT) :: lsum(nsum), 
!!                      integer(IN) :: nsum,
!!                      integer(IN) :: blockID,
!!                      integer(IN) :: idensvar)
!!
!!
!! DESCRIPTION
!!
!!  Computes a block's contributions to some special global sums;
!!  accumulates these contributions into the array lsum.
!!  
!! ARGUMENTS
!!
!!  lsum --     intent(INOUT) array to hold values of special sums
!!         lsum(1) holds a sum of the density weighted by the area of the grid face
!!         lsum(2) holds a lsum(1) * length of x axis
!!         lsum(3) holds a lsum(1) * length of y axis
!!         lsum(4) holds a lsum(1) * length of z axis
!!  nsum --     dimension of the output array (which certainly should be more than 4)
!!  blockID --  current local block index
!!  idensvar -- variable index (from Flash.h) for the grid density variable to be summerd
!!  
!! NOTES
!!  Special version for isolated boundary image mass.
!!
!!***

!!REORDER(4): solnData

subroutine gr_isoSumLocal (lsum, nsum, blockID, idensvar)

!==============================================================================

  use gr_isoMpoleData, ONLY: mpole_geometry, cylfactor, mpole_quadrant, fourpi
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY: Grid_getBlkPtr, Grid_releaseBlkPtr, Grid_getBlkIndexLimits, &
                            Grid_getDeltas, Grid_getBlkBC, Grid_getBlkBoundBox

  implicit none

#include "constants.h"
#include "Flash.h"
  
  integer,intent(IN) :: blockID, idensvar, nsum
  real,intent(INOUT) :: lsum(nsum)
  
  integer,dimension(LOW:HIGH,MDIM) :: faces,blkLimits,blkLimitsGC
  real,dimension(LOW:HIGH,MDIM) :: bnd_box
  real               :: xx, yy, zz
  real, pointer      :: solnData(:,:,:,:)
  integer            :: i, j, k, n
  real               :: dvol, delx, dely, delz, delm

  real               :: delta(MDIM)

  !=====================================================================

  !               Compute dimensions of each zone.
  call Grid_getDeltas(blockID,delta)
  delx = delta(IAXIS)
  if (NDIM >= 2) dely = delta(JAXIS)
  if (NDIM == 3) delz = delta(KAXIS)

  !               Obtain a pointer to the block and its neighbor list.
  call Grid_getBlkPtr(blockID,solnData,CENTER)
  call Grid_getBlkBC(blockID,faces)
  call Grid_getBlkBoundBox(blockID,bnd_box)
  !               Sum contributions from this block.
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  
  select case (mpole_geometry)
     
  case (CARTESIAN)
     if(NDIM==3) then
        if (faces(LOW,IAXIS) /= NOT_BOUNDARY) then
           xx = bnd_box(LOW,IAXIS)
           do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
              zz = bnd_box(LOW,KAXIS) + (k-blkLimits(LOW,KAXIS)+0.5)*delz
              do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
                 yy = bnd_box(LOW,JAXIS) + (j-blkLimits(LOW,JAXIS)+0.5)*dely
                 delm    = solnData(idensvar,blkLimits(LOW,IAXIS)-1,j,k)&
                      * dely * delz
                 lsum(1) = lsum(1) + delm
                 lsum(2) = lsum(2) + delm*xx
                 lsum(3) = lsum(3) + delm*yy
                 lsum(4) = lsum(4) + delm*zz
              enddo
           enddo
        endif
        
        if (faces(HIGH,IAXIS) /= NOT_BOUNDARY) then
           xx = bnd_box(HIGH,IAXIS) 
           do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
              zz = bnd_box(LOW,KAXIS) + (k-blkLimits(LOW,KAXIS)+0.5)*delz
              do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
                 yy = bnd_box(LOW,JAXIS) + (j-blkLimits(LOW,JAXIS)+0.5)*dely
                 delm    = solnData(idensvar,blkLimits(HIGH,IAXIS)+1,j,k) &
                      * dely * delz
                 lsum(1) = lsum(1) + delm
                 lsum(2) = lsum(2) + delm*xx
                 lsum(3) = lsum(3) + delm*yy
                 lsum(4) = lsum(4) + delm*zz
              enddo
           enddo
        endif
        
        if (faces(LOW,JAXIS) /= NOT_BOUNDARY) then
           yy = bnd_box(LOW,JAXIS)
           do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
              zz = bnd_box(LOW,KAXIS) + (k-blkLimits(LOW,KAXIS)+0.5)*delz
              do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
                 xx = bnd_box(LOW,IAXIS) + (i-blkLimits(LOW,IAXIS)+0.5)*delx
                 delm    = solnData(idensvar,i,blkLimits(LOW,JAXIS)-1,k) * delx * delz
                 lsum(1) = lsum(1) + delm
                 lsum(2) = lsum(2) + delm*xx
                 lsum(3) = lsum(3) + delm*yy
                 lsum(4) = lsum(4) + delm*zz
              enddo
           enddo
        endif
        
        if (faces(HIGH,JAXIS) /= NOT_BOUNDARY) then
           yy = bnd_box(HIGH,JAXIS) 
           do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
              zz = bnd_box(LOW,KAXIS) + (k-blkLimits(LOW,KAXIS)+0.5)*delz
              do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
                 xx = bnd_box(LOW,IAXIS) + (i-blkLimits(LOW,IAXIS)+0.5)*delx
                 delm    = solnData(idensvar,i,blkLimits(HIGH,JAXIS)+1,k) * delx * delz
                 lsum(1) = lsum(1) + delm
                 lsum(2) = lsum(2) + delm*xx
                 lsum(3) = lsum(3) + delm*yy
                 lsum(4) = lsum(4) + delm*zz
              enddo
           enddo
        endif
        
        if (faces(LOW,KAXIS) /= NOT_BOUNDARY) then
           zz = bnd_box(LOW,KAXIS)
           do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
              yy = bnd_box(LOW,JAXIS) + (j-blkLimits(LOW,JAXIS)+0.5)*dely
              do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
                 xx = bnd_box(LOW,IAXIS) + (i-blkLimits(LOW,IAXIS)+0.5)*delx
                 delm    = solnData(idensvar,i,j,blkLimits(LOW,KAXIS)-1) * delx * dely
                 lsum(1) = lsum(1) + delm
                 lsum(2) = lsum(2) + delm*xx
                 lsum(3) = lsum(3) + delm*yy
                 lsum(4) = lsum(4) + delm*zz
              enddo
           enddo
        endif
        
        if (faces(HIGH,KAXIS) /= NOT_BOUNDARY) then
           zz = bnd_box(HIGH,KAXIS) 
           do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
              yy = bnd_box(LOW,JAXIS) + (j-blkLimits(LOW,JAXIS)+0.5)*dely
              do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
                 xx = bnd_box(LOW,IAXIS) + (i-blkLimits(LOW,IAXIS)+0.5)*delx
                 delm    = solnData(idensvar,i,j,blkLimits(HIGH,KAXIS)+1) * delx * dely
                 lsum(1) = lsum(1) + delm
                 lsum(2) = lsum(2) + delm*xx
                 lsum(3) = lsum(3) + delm*yy
                 lsum(4) = lsum(4) + delm*zz
              enddo
           enddo
        endif
     else
        call Driver_abortFlash("Isobnd_mpole: cartesian works 3D only")
     end if
     
  case (CYLINDRICAL)
     if(NDIM==2) then
        if (faces(HIGH,IAXIS) /= NOT_BOUNDARY) then
           xx = bnd_box(HIGH,IAXIS) 
           do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
              yy = bnd_box(LOW,JAXIS) + (j-blkLimits(LOW,JAXIS)+0.5)*dely
              delm = solnData(idensvar,blkLimits(HIGH,IAXIS)+1,j,1) * cylfactor * xx * dely
              lsum(1) = lsum(1) + delm
              ! center of mass is on x-axis in 2D cylindrical
              ! center of mass is on y-axis if doing only a quadrant
              if (.not. mpole_quadrant) lsum(3) = lsum(3) + delm*yy
           enddo
        endif
        
        if (faces(LOW,JAXIS) /= NOT_BOUNDARY) then
           yy = bnd_box(LOW,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              xx = bnd_box(LOW,IAXIS) + (i-blkLimits(LOW,IAXIS)+0.5)*delx
              delm = solnData(idensvar,i,blkLimits(LOW,JAXIS)-1,1) * cylfactor * xx * delx
              lsum(1) = lsum(1) + delm
              ! center of mass is on x-axis in 2D cylindrical
              ! center of mass is on y-axis if doing only a quadrant
              if (.not. mpole_quadrant) lsum(3) = lsum(3) + delm*yy
           enddo
        endif
        
        if (faces(HIGH,JAXIS) /= NOT_BOUNDARY) then
           yy = bnd_box(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              xx = bnd_box(LOW,IAXIS) + (i-blkLimits(LOW,IAXIS)+0.5)*delx
              delm = solnData(idensvar,i,blkLimits(LOW,JAXIS)-1,1) * cylfactor * xx * delx
              lsum(1) = lsum(1) + delm
              ! center of mass is on x-axis in 2D cylindrical
              ! center of mass is on y-axis if doing only a quadrant
              if (.not. mpole_quadrant) lsum(3) = lsum(3) + delm*yy
           enddo
        endif
     else
        call Driver_abortFlash("Isobnd_mpole: cylindrical works 2D only")
     end if
     
  case (SPHERICAL)
     if(NDIM==1) then
        if (faces(HIGH,IAXIS) /= NOT_BOUNDARY) then
           xx = bnd_box(HIGH,IAXIS) 
           delm = solnData(idensvar,blkLimits(HIGH,IAXIS)+1,1,1) * fourpi * xx**2
           lsum(1) = lsum(1) + delm
           ! center of mass is at origin in 1D spherical
        endif
     else
        call Driver_abortFlash("Isobnd_mpole: spherical works 1D only")
     end if
     
  end select
  
  call Grid_releaseBlkPtr(blockID,solnData)

  !========================================================================
  
  return
end subroutine gr_isoSumLocal
