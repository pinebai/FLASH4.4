!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserBeams/ed_beamsCheck3DRec
!!
!! NAME
!!
!!  ed_beamsCheck3DRec
!!
!! SYNOPSIS
!!
!!  call ed_beamsCheck3DRec ()
!!
!! DESCRIPTION
!!
!!  Checks the collected beams data for those geometries consisting formally of 3D rectangular
!!  grids (cartesian). All checks which depend on domain grid details should go in here.
!!  Currently it contains the following:
!!
!!   1) Check, if all beam elliptical/rectangular target areas are completely within the domain.
!!   2) Check, if all beam elliptical/rectangular lens areas are completely outside the domain.
!!
!! ARGUMENTS
!!
!! NOTES
!!
!!***

subroutine ed_beamsCheck3DRec ()

  use Driver_interface,         ONLY : Driver_abortFlash

  use EnergyDeposition_data,    ONLY : ed_beams,         &
                                       ed_numberOfBeams, &
                                       ed_xminDomain,    &
                                       ed_xmaxDomain,    &
                                       ed_yminDomain,    &
                                       ed_ymaxDomain,    &
                                       ed_zminDomain,    &
                                       ed_zmaxDomain
  
  implicit none

#include "Flash.h"
#include "EnergyDeposition.h"
#include "constants.h"

  character (len = BEAM_STRING_LENGTH) :: gridType

  logical :: inDomain
  logical :: outOfDomain

  integer :: beam

  real    :: c1, c2, c3, c4
  real    :: cxmin, cxmax, cymin, cymax, czmin, czmax
  real    :: ex, ey, ez
  real    :: lensX, lensY, lensZ
  real    :: s1, s2
  real    :: targetX, targetY, targetZ
  real    :: u1x, u1y, u1z
  real    :: u2x, u2y, u2z
!
!
!     ...Loop over all beams.
!
!
  do beam = 1, ed_numberOfBeams

     lensX    = ed_beams (beam) % lensX
     lensY    = ed_beams (beam) % lensY
     lensZ    = ed_beams (beam) % lensZ
     gridType = ed_beams (beam) % gridType
     u1x      = ed_beams (beam) % semiAxisUnitMajorX
     u1y      = ed_beams (beam) % semiAxisUnitMajorY
     u1z      = ed_beams (beam) % semiAxisUnitMajorZ
     u2x      = ed_beams (beam) % semiAxisUnitMinorX
     u2y      = ed_beams (beam) % semiAxisUnitMinorY
     u2z      = ed_beams (beam) % semiAxisUnitMinorZ
     targetX  = ed_beams (beam) % targetX
     targetY  = ed_beams (beam) % targetY
     targetZ  = ed_beams (beam) % targetZ

     if (gridType /= 'rectangular2D') then
!
!
!     ...The elliptical crossection case:
!        -------------------------------
!
!        The implicit form of the elliptical target boundary curve is:
!
!              e = s1 * u1 * cos (t)  +  s2 * u2 * sin (t)
!
!        where u1,u2 stand for the two semiaxis unit vectors in the local target coordinate
!        system, s1,s2 are the two semiaxes lengths, 'e' is a vector on the elliptical boundary
!        and 't' is the implicit parameter ranging from 0 to 2pi. Differentiating this
!        equation with respect to 't' we obtain:
!
!                     t = arctan ( [s2 * u2] / [s1 * u1] )
!
!        which, when simplifying, leads to the simple minimax equation:
!
!                  e (min,max) = +/- sqrt ([s1 *u1]^2 + [s2 *u2]^2)
!
!        or, in component form:
!
!                 ex (min,max) = +/- sqrt ([s1 *u1x]^2 + [s2 *u2x]^2)
!                 ey (min,max) = +/- sqrt ([s1 *u1y]^2 + [s2 *u2y]^2)
!                 ez (min,max) = +/- sqrt ([s1 *u1z]^2 + [s2 *u2z]^2)
!
!        Since the domain boundaries are in terms of the global coordinate system, we need
!        to convert the elliptical curve points from the local to the global coordinate
!        system via:
!
!                                E = e + T
!
!        where T is the position of the target center. The results of E will then
!        be tested against the domain boundaries.
!
!
         s1 = ed_beams (beam) % targetSemiAxisMajor
         s2 = ed_beams (beam) % targetSemiAxisMinor

         ex = sqrt ( (s1 * u1x) ** 2 + (s2 * u2x) ** 2 )
         ey = sqrt ( (s1 * u1y) ** 2 + (s2 * u2y) ** 2 )
         ez = sqrt ( (s1 * u1z) ** 2 + (s2 * u2z) ** 2 )
!
!
!     ...Check, if beam elliptical target area is incident completely within domain.
!
!
         inDomain =      (targetX - ex >= ed_xminDomain) &
                   .and. (targetY - ey >= ed_yminDomain) &
                   .and. (targetZ - ez >= ed_zminDomain) &
                   .and. (targetX + ex <= ed_xmaxDomain) &
                   .and. (targetY + ey <= ed_ymaxDomain) &
                   .and. (targetZ + ez <= ed_zmaxDomain)

         if (.not.inDomain) then
             call Driver_abortFlash ("ed_beamsCheck3DRec: Beam target (partially) outside of domain!")
         end if
!
!
!     ...Check, if beam elliptical lens area is completely outside of the domain.
!
!
         s1 = ed_beams (beam) % lensSemiAxisMajor
         s2 = ed_beams (beam) % lensSemiAxisMinor

         ex = sqrt ( (s1 * u1x) ** 2 + (s2 * u2x) ** 2 )
         ey = sqrt ( (s1 * u1y) ** 2 + (s2 * u2y) ** 2 )
         ez = sqrt ( (s1 * u1z) ** 2 + (s2 * u2z) ** 2 )

         outOfDomain =     (lensX - ex > ed_xmaxDomain) &
                      .or. (lensY - ey > ed_ymaxDomain) &
                      .or. (lensZ - ez > ed_zmaxDomain) &
                      .or. (lensX + ex < ed_xminDomain) &
                      .or. (lensY + ey < ed_yminDomain) &
                      .or. (lensZ + ez < ed_zminDomain)

         if (.not.outOfDomain) then
             call Driver_abortFlash ("ed_beamsCheck3DRec: Beam lens (partially) inside the domain!")
         end if

     else
!
!
!     ...The rectangular crossection case:
!        --------------------------------
!
!        Since the target and lens areas are rectangular, the extremum values in each direction
!        correspond to one of the four corners. Using the target/lens center T/L global X,Y,Z
!        coordinates and the local target/lens x,y,z coordinates, each corner C has the following
!        global coordinate:
!
!                 C(X,Y,Z) = T/L(X,Y,Z) +/- s1 * u1(x,y,z) +/- s2 * u2(x,y,z)]
!
!        where u1,u2 stand for the two semiaxis unit vectors in the local target/lens coordinate
!        system and s1,s2 are the two corresponding semiaxes lengths:
!
!
!                                    C\
!                                    | \
!                                    |  \
!                                    |   \
!                                    | u1|\
!                                    |   | \        T/L = target/lens center
!                                    |   |  \
!                                    |   |   C        C = target/lens corners
!                                    |  T/L  |
!                                    C    \  |    u1,u2 = target/lens unit vectors
!                                     \    \ |
!                                      \  u2\|
!                                       \    |
!                                        \   |
!                                         \  |
!                                          \ |
!                                           C|
!
!
!        For each X,Y,Z component of the corners we need to find the maximum and the minimum
!        between the four corner values and compare it to the domain extensions.
!        First for the target.
!
!
         s1 = ed_beams (beam) % targetSemiAxisMajor
         s2 = ed_beams (beam) % targetSemiAxisMinor

         c1 = targetX + s1 * u1x + s2 * u2x
         c2 = targetX + s1 * u1x - s2 * u2x
         c3 = targetX - s1 * u1x + s2 * u2x
         c4 = targetX - s1 * u1x - s2 * u2x

         cxmin = min (c1,c2,c3,c4)
         cxmax = max (c1,c2,c3,c4)

         c1 = targetY + s1 * u1y + s2 * u2y
         c2 = targetY + s1 * u1y - s2 * u2y
         c3 = targetY - s1 * u1y + s2 * u2y
         c4 = targetY - s1 * u1y - s2 * u2y

         cymin = min (c1,c2,c3,c4)
         cymax = max (c1,c2,c3,c4)

         c1 = targetZ + s1 * u1z + s2 * u2z
         c2 = targetZ + s1 * u1z - s2 * u2z
         c3 = targetZ - s1 * u1z + s2 * u2z
         c4 = targetZ - s1 * u1z - s2 * u2z

         czmin = min (c1,c2,c3,c4)
         czmax = max (c1,c2,c3,c4)
!
!
!     ...Check, if beam rectangular target area is incident completely within domain.
!
!
         inDomain =      (cxmin >= ed_xminDomain) &
                   .and. (cymin >= ed_yminDomain) &
                   .and. (czmin >= ed_zminDomain) &
                   .and. (cxmax <= ed_xmaxDomain) &
                   .and. (cymax <= ed_ymaxDomain) &
                   .and. (czmax <= ed_zmaxDomain)

         if (.not.inDomain) then
             call Driver_abortFlash ("ed_beamsCheck3DRec: Beam target (partially) outside of domain!")
         end if
!
!
!     ...Check, if beam rectangular lens area is completely outside of the domain.
!
!
         s1 = ed_beams (beam) % lensSemiAxisMajor
         s2 = ed_beams (beam) % lensSemiAxisMinor

         c1 = lensX + s1 * u1x + s2 * u2x
         c2 = lensX + s1 * u1x - s2 * u2x
         c3 = lensX - s1 * u1x + s2 * u2x
         c4 = lensX - s1 * u1x - s2 * u2x

         cxmin = min (c1,c2,c3,c4)
         cxmax = max (c1,c2,c3,c4)

         c1 = lensY + s1 * u1y + s2 * u2y
         c2 = lensY + s1 * u1y - s2 * u2y
         c3 = lensY - s1 * u1y + s2 * u2y
         c4 = lensY - s1 * u1y - s2 * u2y

         cymin = min (c1,c2,c3,c4)
         cymax = max (c1,c2,c3,c4)

         c1 = lensZ + s1 * u1z + s2 * u2z
         c2 = lensZ + s1 * u1z - s2 * u2z
         c3 = lensZ - s1 * u1z + s2 * u2z
         c4 = lensZ - s1 * u1z - s2 * u2z

         czmin = min (c1,c2,c3,c4)
         czmax = max (c1,c2,c3,c4)

         outOfDomain =     (cxmin > ed_xmaxDomain) &
                      .or. (cymin > ed_ymaxDomain) &
                      .or. (czmin > ed_zmaxDomain) &
                      .or. (cxmax < ed_xminDomain) &
                      .or. (cymax < ed_yminDomain) &
                      .or. (czmax < ed_zminDomain)

         if (.not.outOfDomain) then
             call Driver_abortFlash ("ed_beamsCheck3DRec: Beam lens (partially) inside the domain!")
         end if

     end if

  end do
!
!
!     ...Ready!
!
!
  return
end subroutine ed_beamsCheck3DRec
