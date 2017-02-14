!!****if* source/Grid/GridMain/paramesh/gr_markEllipsoid
!!  
!! NAME 
!!  gr_markEllipsoid 
!!  
!! SYNOPSIS 
!!  gr_markEllipsoid(real (in) :: ic, 
!!                             real (in) :: jc, 
!!                             real (in) :: kc, 
!!                             real (in) :: a1, 
!!                             real (in) :: a2, 
!!                             real (in) :: a3, 
!!                             integer (in) :: lref) 
!!  
!! PURPOSE 
!!  Refine all blocks containing points on an ellipsoidal surface centered on
!!  (ic,jc,kc) with semimajor axes (a1,a2,a3).  Either blocks are brought up to
!!  a specific level of refinement or each block is refined once.  
!!  
!! ARGUMENTS 
!!  ic - i coordinate of the ellipsoid center
!!  jc - j coordinate of the ellipsoid center
!!  kc - k coordinate of the ellipsoid center
!!               (Coordinates for nonexistent dimensions are ignored.)
!!  a1 - Semimajor axis of the ellipsoid
!!  a2 - Semimajor axis of the ellipsoid
!!  a3 - Semimajor axis of the ellipsoid
!!               (Axes for nonexistent dimensions are ignored.)
!!  lref-        If > 0, bring all qualifying blocks to this level of refinement.
!!
!!               If <= 0, refine qualifying blocks once.
!!  
!! NOTES
!! 
!!  This routine has not yet been tested and should be used only as a guideline for
!!  a user's implementation.
!!  
!!***

subroutine gr_markEllipsoid(ic, jc, kc, a1, a2, a3, lref)

!-------------------------------------------------------------------------------
  use tree, ONLY : refine, derefine, lrefine, coord, bsize, lnblocks, nodetype
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_data, ONLY : gr_geometry

  implicit none
#include "constants.h"
#include "Flash.h"
! Arguments

  real, intent(IN)      :: ic, jc, kc, a1, a2, a3
  integer, intent(IN)   :: lref

! Local data

  real, dimension(MDIM) :: blockCenter, blockSize
  real                  :: bxl, bxr, byl, byr, bzl, bzr
  real                  :: mindist2, xmindist2, ymindist2, zmindist2
  real                  :: maxdist2, xmaxdist2, ymaxdist2, zmaxdist2
  real                  :: a1inv2, a2inv2, a3inv2
  integer               :: b



!-------------------------------------------------------------------------------
  if((gr_geometry /= CARTESIAN).and.(gr_geometry /= CYLINDRICAL)) &
       call Driver_abortFlash("markrefineEllipsoid : wrong geometry")
  if((gr_geometry == CYLINDRICAL).and.(NDIM==3))&
       call Driver_abortFlash("markrefineEllipsoid : 3d not supported in cylindrical")
! Compute squared inverses of semi-major axes

  a1inv2 = 1./a1**2
  a2inv2 = 1./a2**2
  a3inv2 = 1./a3**2

! Loop over leaf-node blocks


  do b = 1, lnblocks
     if (nodetype(b) == 1) then

! Get bounding box for this block (relative to (ic,jc,kc))

        blockCenter = coord(:,b)
        blockSize  = 0.5 * bsize(:,b)
        
        bxl = blockCenter(1) - blockSize(1) - ic
        bxr = blockCenter(1) + blockSize(1) - ic
        if (NDIM > 1) then
           byl = blockCenter(2) - blockSize(2) - jc
           byr = blockCenter(2) + blockSize(2) - jc
        else
           byl = 0.
           byr = 0.
        endif
        if (NDIM == 3) then
           bzl = blockCenter(3) - blockSize(3) - kc
           bzr = blockCenter(3) + blockSize(3) - kc
        else
           bzl = 0.
           bzr = 0.
        endif
        
! Find minimum and maximum distance from (ic,jc,kc) for each dimension.
! For each coordinate, if both "left" and "right" distances have the same sign,
! then the smaller magnitude is the minimum.  Otherwise (ic,jc,kc) is
! contained within the interval for that dimension, so the minimum is 0.
! The maximum distance is always the larger of the two magnitudes.
! Nonexistent dimensions have had all distances set to zero, so they are
! ignored.

! For an ellipsoidal surface, scale all distances by the semimajor axes.

        if (bxl*bxr > 0.) then
           xmindist2 = min( bxl**2, bxr**2 ) * a1inv2
        else
           xmindist2 = 0.
        endif
        xmaxdist2 = max( bxl**2, bxr**2 ) * a1inv2
        
        if (byl*byr > 0.) then
           ymindist2 = min( byl**2, byr**2 ) * a2inv2
        else
           ymindist2 = 0.
        endif
        ymaxdist2 = max( byl**2, byr**2 ) * a2inv2
        
        if (bzl*bzr > 0.) then
           zmindist2 = min( bzl**2, bzr**2 ) * a3inv2
        else
           zmindist2 = 0.
        endif
        zmaxdist2 = max( bzl**2, bzr**2 ) * a3inv2
        
        ! Now compute the minimum and maximum scaled distances to (ic,jc,kc) and
        ! compare them to 1.  If unity is within the interval defined by the
        ! two scaled distances, then the block contains points on the surface of the
        ! ellipsoid and is marked for refinement.
        
        mindist2 = xmindist2 + ymindist2 + zmindist2
        maxdist2 = xmaxdist2 + ymaxdist2 + zmaxdist2
        ! Currently assumes Cartesian
        ! or 2D axisymmetric (r-z)
        ! or 1D spherical (r)
        if ((maxdist2 >= 1.) .and. (mindist2 <= 1.)) then
           
           if (lrefine(b) < lref ) then
              refine(b)   = .true.
              derefine(b) = .false.
           else if (lrefine(b) == lref) then
              derefine(b) = .false.
           else if (lref <= 0) then
              refine(b) = .true.
           endif
           
        endif
        
        ! End of leaf-node block loop
        
     endif
  end do
  
  !-------------------------------------------------------------------------------
  
  return
end subroutine gr_markEllipsoid
