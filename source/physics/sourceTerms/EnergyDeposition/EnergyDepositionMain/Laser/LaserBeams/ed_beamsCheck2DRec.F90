!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserBeams/ed_beamsCheck2DRec
!!
!! NAME
!!
!!  ed_beamsCheck2DRec
!!
!! SYNOPSIS
!!
!!  call ed_beamsCheck2DRec ()
!!
!! DESCRIPTION
!!
!!  Checks the collected beams data for those geometries consisting formally of 2D rectangular
!!  grids (cartesian + cylindrical). All checks which depend on domain grid details should go
!!  in here. Currently it contains the following:
!!
!!         1) Check, if all beam line target areas are completely within the domain.
!!         2) Check, if all beam line lens areas are completely outside the domain.
!!
!! ARGUMENTS
!!
!! NOTES
!!
!!***

subroutine ed_beamsCheck2DRec ()

  use Driver_interface,         ONLY : Driver_abortFlash

  use EnergyDeposition_data,    ONLY : ed_beams,         &
                                       ed_gridGeometry,  &
                                       ed_numberOfBeams, &
                                       ed_xminDomain,    &
                                       ed_xmaxDomain,    &
                                       ed_yminDomain,    &
                                       ed_ymaxDomain
  
  implicit none

#include "Flash.h"
#include "EnergyDeposition.h"
#include "constants.h"

  logical :: inDomain
  logical :: outOfDomain

  integer :: beam

  real    :: lensX, lensY
  real    :: Lx, Ly
  real    :: s
  real    :: targetX, targetY
  real    :: ux, uy
  real    :: xmaxDomain, xminDomain
  real    :: ymaxDomain, yminDomain
!
!
!     ...In case of 2D cylindrical geometry, treat the 2D cylindrical domain as the fully
!        3-dimensional rotational domain. This allows placing lenses and targets in the
!        negative x-coordinates area.
!
!
  xmaxDomain = ed_xmaxDomain
  yminDomain = ed_yminDomain
  ymaxDomain = ed_ymaxDomain

  if (ed_gridGeometry == GRID_2DCYLINDRICAL) then
      xminDomain = - ed_xmaxDomain
  else
      xminDomain =   ed_xminDomain
  end if
!
!
!     ...Loop over all beams.
!
!
  do beam = 1, ed_numberOfBeams

     ux      = ed_beams (beam) % semiAxisUnitMajorX
     uy      = ed_beams (beam) % semiAxisUnitMajorY
     s       = ed_beams (beam) % targetSemiAxisMajor
     targetX = ed_beams (beam) % targetX
     targetY = ed_beams (beam) % targetY
!
!
!     ...Any point on a line along the target unit vectors is given by:
!
!                               L = t * u
!
!        where 'u' stand for the line unit vector in the local target coordinate
!        system, 't' is a parameter from - infinity to + infinity and 'L' is a vector
!        on the infinite line containing the target line. The maximum values that
!        Lx and Ly can attend for the target are simply:
!
!                          L (min,max) = +/- s * u
!
!        where 's' is half the target line length. In component form we have:
!
!                         Lx (min,max) = +/- s * ux
!                         Ly (min,max) = +/- s * uy
!
!        Since the domain boundaries are in terms of the global coordinate system, we need
!        to convert the target line points from the local to the global coordinate
!        system via forming 'L + T', where T is the position of the target center.
!        The results of 'L + T' will then be tested against the domain boundaries.
!
!
     Lx = s * abs (ux)
     Ly = s * abs (uy)
!
!
!     ...Check, if beam target area is incident completely within domain.
!
!
     inDomain =      (targetX - Lx >= xminDomain) &
               .and. (targetY - Ly >= yminDomain) &
               .and. (targetX + Lx <= xmaxDomain) &
               .and. (targetY + Ly <= ymaxDomain)

     if (.not.inDomain) then
         call Driver_abortFlash ("ed_beamsCheck2DRec: Beam target (partially) outside of domain!")
     end if
!
!
!     ...Check, if beam lens area is completely outside of the domain.
!
!
     s     = ed_beams (beam) % lensSemiAxisMajor
     lensX = ed_beams (beam) % lensX
     lensY = ed_beams (beam) % lensY

     Lx = s * abs (ux)
     Ly = s * abs (uy)

     outOfDomain =     (lensX - Lx > xmaxDomain) &
                  .or. (lensY - Ly > ymaxDomain) &
                  .or. (lensX + Lx < xminDomain) &
                  .or. (lensY + Ly < yminDomain)

     if (.not.outOfDomain) then
         call Driver_abortFlash ("ed_beamsCheck2DRec: Beam lens (partially) inside the domain!")
     end if

  end do
!
!
!     ...Ready!
!
!
  return
end subroutine ed_beamsCheck2DRec
