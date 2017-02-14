!!****if* source/Grid/GridMain/paramesh/gr_markWithRadius
!!  
!! NAME 
!!  gr_markWithRadius 
!!  
!! SYNOPSIS 
!!  gr_markWithRadius(real(IN) :: ic, 
!!                       real(IN) :: jc, 
!!                       real(IN) :: kc, 
!!                       real(IN) :: radius, 
!!                       integer(IN) :: lref) 
!!  
!! PURPOSE 
!!  Refine all blocks containing points at a given distance from a given point
!!  (ic,jc,kc).  Either blocks are brought up to a specific level of refinement
!!  or each block is refined once.  
!!  
!! ARGUMENTS 
!!  ic, jc, kc:  Center of the interval/circle/sphere
!!               (Coordinates for nonexistent dimensions are ignored.)
!!  radius:      Radius of the region 
!!   lref -       If > 0, bring all qualifying blocks to this level of refinement.
!!
!!               If <= 0, refine qualifying blocks once.
!!
!! NOTES
!! 
!!  This routine has not yet been tested and should be used only as a guideline for
!!  a user's implementation.
!!  
!!  
!!***

subroutine gr_markWithRadius(ic, jc, kc, radius, lref)

!-------------------------------------------------------------------------------
  use tree, ONLY : refine, derefine, lrefine, nodetype, lnblocks, coord, bsize
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_data, ONLY : gr_geometry
  implicit none
#include "constants.h"
#include "Flash.h"

  real, intent(IN)      :: ic, jc, kc, radius
  integer, intent(IN)   :: lref

! Local data

  real, dimension(MDIM) :: blockCenter, blockSize
  real                  :: bxl, bxr, byl, byr, bzl, bzr
  real                  :: mindist2, xmindist2, ymindist2, zmindist2
  real                  :: maxdist2, xmaxdist2, ymaxdist2, zmaxdist2
  integer               :: b

!------------------------------------------------------------------------------

  if((gr_geometry == CARTESIAN).or.(gr_geometry == CYLINDRICAL)) then
     do b = 1, lnblocks
        if (nodetype(b) == LEAF) then

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
           if ((NDIM == 3).and.(gr_geometry==CYLINDRICAL)) then
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

           if (bxl*bxr > 0.) then
              xmindist2 = min( bxl**2, bxr**2 )
           else
              xmindist2 = 0.
           endif
           xmaxdist2 = max( bxl**2, bxr**2 )
           
           if (byl*byr > 0.) then
              ymindist2 = min( byl**2, byr**2 )
           else
              ymindist2 = 0.
           endif
           ymaxdist2 = max( byl**2, byr**2 )
           
           if (bzl*bzr > 0.) then
              zmindist2 = min( bzl**2, bzr**2 )
           else
              zmindist2 = 0.
           endif
           zmaxdist2 = max( bzl**2, bzr**2 )
           
! Now compute the minimum and maximum distances to (ic,jc,kc) and compare them
! to the specified radius.  If the radius is within the interval defined by the
! two distances, then the block contains points on the boundary of the
! interval/circle/sphere and is Grid_marked for refinement.

           mindist2 = xmindist2 + ymindist2 + zmindist2
           maxdist2 = xmaxdist2 + ymaxdist2 + zmaxdist2
           ! Currently assumes Cartesian
           ! or 2D axisymmetric (r-z)
           ! or 1D spherical (r)
           if ((maxdist2 >= radius**2) .and. (mindist2 <= radius**2)) then

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
  elseif((gr_geometry==POLAR).or.(gr_geometry==SPHERICAL)) then
     do b = 1, lnblocks
        if (nodetype(b) == LEAF) then

! Get bounding box for this block (relative to (ic,jc,kc))

           blockCenter = coord(:,b)
           blockSize  = 0.5 * bsize(:,b)
           
           bxl = blockCenter(1) - blockSize(1) - ic
           bxr = blockCenter(1) + blockSize(1) - ic
           if (bxl*bxr > 0.) then
              mindist2 = min( bxl, bxr)
           else
              mindist2 = 0.
           endif
           maxdist2 = max( bxl, bxr )
           if ((maxdist2 >= radius) .and. (mindist2 <= radius)) then

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
  else
     call Driver_abortFlash("gr_markWithRadius : wrong geometry")
  end if
  !-------------------------------------------------------------------------------
  
  return
end subroutine gr_markWithRadius
