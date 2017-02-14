!!****if* source/Grid/GridSolvers/Multipole/gr_mpoleLocalSum
!!
!! NAME
!!
!!  gr_mpoleLocalSum
!!
!! SYNOPSIS
!!
!!  call gr_mpoleLocalSum(integer(in) :: blockID,
!!                        integer(in) :: nsum,
!!                        integer(in) :: idensvar,
!!                        real (out)  :: lsum(nsum))
!!
!! DESCRIPTION
!!
!!   Accumulate values of the density and density-weighted position
!!   from the interior of a given block.  The results are stored in
!!   a given vector (lsum), which is assumed to have been
!!   initialized by the calling routine.
!!
!! ARGUMENTS
!!
!!   blockID - the block
!!   nsum    - size of lsum array
!!   idensvar- the index for density in mesh data structure
!!   lsum    - array that returns accumulated values
!!
!!***

!!REORDER(4): solnData

subroutine gr_mpoleLocalSum (blockID,nsum,idensvar, lsum)
  
  !==================================================================
  
  use gr_mpoleData, ONLY : G_3DCARTESIAN,G_1DSPHERICAL,G_2DCYLINDRICAL,G_3DAXISYMMETRIC,&
                         quadrant, twopi, fourpi, mpole_geometry,&
                         octant
  use Grid_interface, ONLY : Grid_getBlkPtr,Grid_releaseBlkPtr,&
                             Grid_getBlkBoundBox,Grid_getDeltas,&
                             Grid_getBlkIndexLimits

  implicit none
  
#include "constants.h"
#include "Flash.h"

  integer,intent(IN) :: blockID, idensvar, nsum
  real,intent(INOUT)  :: lsum(nsum)
  
  real,dimension(MDIM) :: delta
  real,dimension(LOW:HIGH,MDIM) :: bndBox
  integer,dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC

  real               :: xx, yy, zz
  real, pointer, dimension(:,:,:,:)      :: solnData
  integer            :: i, j, k, imax, jmax, kmax, imin, jmin, kmin
  real               :: dvol,  delm
  
  !=====================================================================
  
  call Grid_getBlkBoundBox(blockID, bndBox)

  ! Compute dimensions of each zone.
  call Grid_getDeltas(blockID,delta)
  call Grid_getBlkPtr(blockID, solnData)

  
  !               Sum contributions from this block.
  
  yy = 0.
  zz = 0.
  
  call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)


  kmax = blkLimits(HIGH,KAXIS)
  kmin = blkLimits(LOW,KAXIS)  
  jmax = blkLimits(HIGH,JAXIS)
  jmin = blkLimits(LOW,JAXIS)
  imax = blkLimits(HIGH,IAXIS)
  imin = blkLimits(LOW,IAXIS)



  do k = kmin, kmax
     if (NDIM == 3) zz = bndBox(LOW,KAXIS) + (k-kmin+0.5)*delta(KAXIS)
     do j = jmin, jmax
        if (NDIM >= 2) yy = bndBox(LOW,JAXIS) + (j-jmin+0.5)*delta(JAXIS)
        do i = imin, imax
           xx = bndBox(LOW,IAXIS) + (i-imin+0.5)*delta(IAXIS)
           
           select case (mpole_geometry)
              
           case (G_3DCARTESIAN)
              dvol = delta(IAXIS) * delta(JAXIS) * delta(KAXIS)
              
              if (octant) then
                 
                 ! if we are doing an octant, force the center
                 ! of mass to be at the origin (due to symmetry)
                 xx = 0.0
                 yy = 0.0
                 zz = 0.0
                 
                 ! we want the mass computed to reflect the total
                 ! mass of the star
                 dvol = 8.0 * dvol
              endif

           case (G_3DAXISYMMETRIC)
              dvol = delta(IAXIS) * delta(JAXIS) * delta(KAXIS)
              xx = 0.            ! ctr of mass stays on z-axis in axisymmetry
              yy = 0.
              
           case (G_2DCYLINDRICAL)
              dvol = twopi * xx * delta(IAXIS) * delta(JAXIS)
              xx = 0.            ! ctr of mass is on x-axis in 2D cylindrical
              if (quadrant) then
                 yy = 0.          ! ctr of mass is on y-axis if doing a quadrant
                 dvol = 2. * dvol ! account for symmetry-suppressed quadrant
              endif
              
           case (G_1DSPHERICAL)
              dvol = fourpi * xx**2 * delta(IAXIS)
              xx = 0.     ! center of mass is at origin in 1D spherical
              
           end select
           
           delm    = solnData(idensvar,i,j,k)*dvol
           lsum(1) = lsum(1) + delm
           lsum(2) = lsum(2) + delm*xx
           lsum(3) = lsum(3) + delm*yy
           lsum(4) = lsum(4) + delm*zz
           
        enddo
     enddo
  enddo
  
  !==========================================================================
  call Grid_releaseBlkPtr(blockID, solnData)
  
  return
end subroutine gr_mpoleLocalSum
