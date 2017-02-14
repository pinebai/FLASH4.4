!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/oldOneTimeStep/ProtonDetection/pi_recordProtonOnScreen
!!
!! NAME
!!
!!  pi_recordProtonOnScreen
!!
!! SYNOPSIS
!!
!!  call pi_recordProtonOnScreen (real    (in)            :: px,
!!                                real    (in)            :: py,
!!                                real    (in)            :: pz,
!!                                real    (in)            :: vx,
!!                                real    (in)            :: vy,
!!                                real    (in)            :: vz,
!!                                real    (in)            :: Jv,
!!                                real    (in)            :: Kx,
!!                                real    (in)            :: Ky,
!!                                real    (in)            :: Kz,
!!                                integer (in)            :: detector,
!!                                logical (out), optional :: onScreen
!!                                real    (out), optional :: sx,
!!                                real    (out), optional :: sy,
!!                                real    (out), optional :: sz)
!!
!! DESCRIPTION
!!
!!  Records the proton on the specified detector screen. The proton has an initial
!!  position and direction (velocity). The routine records the (x,y) location on
!!  the screen in terms of local screen coordinates as well as certain diagnostic
!!  variables (currently Jv,Kx,Ky,Kz). The (x,y) screen coordinates are rescaled to
!!  the square screen side length, such that all protons within the screen have screen
!!  coordinate paris (x,y) within the range [0,1]. Protons falling outside the screen
!!  are also saved for eventual further diagnostic studies. Protons which miss the
!!  screen plane entirely are not recorded.
!!
!!  Procedure:
!!
!!                                    |\
!!                                    | \
!!                                    |  \                   C = detector center
!!                                    |   \
!!                                    | uY|\                 n = normal unit vector
!!                                    |   | \
!!                                    |   |  \         uX,uY,n = form a right hand rule
!!                      detector      |   |   |                  orthogonal vector system:
!!                                    |   C---|---> n                  n = uY x uX
!!                       screen       \  / \  |
!!                                     \/   \ |              P = position of proton
!!                                     /\  uX\|
!!                                    /  \    |              v = velocity vector of proton
!!                                   /    \   |
!!                                  /      \  |          (x,y) = local screen coordinates
!!                                 /        \ |                  (basis = uX,uY) where proton
!!                                /          \|                   will hit the screen plane
!!                               /
!!                              /                   (sx,sy,sz) = global screen coordinates,
!!                             /                                 calculated only, if proton hits
!!                            /                                  screen plane. (0,0,0) otherwise.
!!                           /
!!                          P---> v ----------- (x,y)
!!                                         (sx,sy,sz)
!!
!!  From geometry:
!!
!!                    i) [(C - P) dot n] / [v dot n] > 0  ->  proton hits screen plane
!!
!!                   ii) x = uY dot [v cross (C - P)] / [v dot n]
!!
!!                  iii) y = uX dot [(C - P) cross v] / [v dot n]
!!
!!                   iv) if [v dot n] = 0  ->  parallel to screen plane, no hit
!!
!!
!!  If a pinhole was attached (non-zero pinhole radius) to the detector screen, the same
!!  equations can be used to check if the proton goes through the pinhole. The current
!!  definition of the pinhole implies that its center H is located opposite to the normal
!!  unit vector n and its circular opening is co-planar with the plane of the detector
!!  screen. Hence the unit vector pair (uX,uY) from the detector screen can be used to
!!  characterize the pinhole circular plane as well. The local pinhole circular plane
!!  coordinates (hx,hy) where the proton will hit the pinhole plane, are calculated in the
!!  same way as for the detector screen, using equations ii) and iii), but with C replaced
!!  by H:
!!  
!!                   ii) hx = uY dot [v cross (H - P)] / [v dot n]
!!
!!                  iii) hy = uX dot [(H - P) cross v] / [v dot n]
!!
!!  From these values the square of the distance hx^2 + hy^2 is formed and checked against
!!  the square of the given pinhole radius, to see if the proton made it through the hole.
!!
!!
!! ARGUMENTS
!!
!!  px       : position global x-coordinate of the proton
!!  py       : position global y-coordinate of the proton
!!  pz       : position global z-coordinate of the proton
!!  vx       : velocity x-coordinate component of the proton
!!  vy       : velocity y-coordinate component of the proton
!!  vz       : velocity z-coordinate component of the proton
!!  Jv       : diagnostic variable
!!  Kx       : diagnostic variable
!!  Ky       : diagnostic variable
!!  Kz       : diagnostic variable
!!  detector : detector screen number on which to record the proton
!!  onScreen : if true, the proton hits the screen
!!  sx       : screen global x-coordinate of the proton (if screen plane is hit, 0 otherwise)
!!  sy       : screen global y-coordinate of the proton (if screen plane is hit, 0 otherwise)
!!  sz       : screen global z-coordinate of the proton (if screen plane is hit, 0 otherwise)
!!
!! NOTES
!!
!!***

subroutine pi_recordProtonOnScreen (px, py, pz,                &
                                    vx, vy, vz,                &
                                    Jv, Kx, Ky, Kz,            &
                                    detector,                  &
                                                    onScreen,  &
                                                    sx, sy, sz )

  use ProtonImaging_data,  ONLY : pi_detectors,              &
                                  pi_recordOffscreenProtons, &
                                  pi_screenProtonCount,      &
                                  pi_screenProtons

  use Driver_interface,    ONLY : Driver_abortFlash

  implicit none

#include "constants.h"
#include "Flash.h"
#include "ProtonImaging.h"

  real,    intent (in)            :: px, py, pz
  real,    intent (in)            :: vx, vy, vz
  real,    intent (in)            :: Jv, Kx, Ky, Kz
  integer, intent (in)            :: detector
  logical, intent (out), optional :: onScreen
  real,    intent (out), optional :: sx, sy, sz

  logical :: hitScreen
  logical :: recordProton

  real    :: cpdotn
  real    :: cpvx, cpvy, cpvz
  real    :: cpx, cpy, cpz
  real    :: cx, cy, cz
  real    :: hpvx, hpvy, hpvz
  real    :: hpx, hpy, hpz
  real    :: hx, hy, hz
  real    :: nx, ny, nz
  real    :: pinholeRadius
  real    :: sHalf, sInv
  real    :: uXx, uXy, uXz
  real    :: uYx, uYy, uYz
  real    :: vdotn, vdotnInv
  real    :: x, y
  real    :: x01, y01
!
!
!     ...Load the detector specifics.
!
!
  nx  = pi_detectors (detector) % normalX
  ny  = pi_detectors (detector) % normalY
  nz  = pi_detectors (detector) % normalZ

  cx  = pi_detectors (detector) % centerX
  cy  = pi_detectors (detector) % centerY
  cz  = pi_detectors (detector) % centerZ

  uXx = pi_detectors (detector) % axisXunitX
  uXy = pi_detectors (detector) % axisXunitY
  uXz = pi_detectors (detector) % axisXunitZ

  uYx = pi_detectors (detector) % axisYunitX
  uYy = pi_detectors (detector) % axisYunitY
  uYz = pi_detectors (detector) % axisYunitZ

  sHalf = pi_detectors (detector) % sideLengthHalf
  sInv  = pi_detectors (detector) % sideLengthInv
!
!
!     ...Form the scalar product of v with n and check, if the proton travels parallel to
!        the screen (a miss). If not, calculate its inverse.
!
!
  vdotn = vx * nx + vy * ny + vz * nz

  if (vdotn == 0.0) then
      if (present (sx)) sx = 0.0
      if (present (sy)) sy = 0.0
      if (present (sz)) sz = 0.0
      if (present (onScreen)) onScreen = .false.
      return
  end if

  vdotnInv = 1.0 / vdotn
!
!
!     ...Form the scalar product of (C-P) with n and check, if the proton travels away from
!        the screen (a miss).
!
!

  cpx = cx - px
  cpy = cy - py
  cpz = cz - pz

  cpdotn = cpx * nx + cpy * ny + cpz * nz

  if (cpdotn * vdotnInv <= 0.0) then
      if (present (sx)) sx = 0.0
      if (present (sy)) sy = 0.0
      if (present (sz)) sz = 0.0
      if (present (onScreen)) onScreen = .false.
      return
  end if
!
!
!     ...We are now sure the proton hits the screen plane. Form the cross product of (c-p)
!        with v and calculate the (x,y) screen coordinates from the center of the screen.
!
!
  cpvx = cpy * vz - cpz * vy
  cpvy = cpz * vx - cpx * vz
  cpvz = cpx * vy - cpy * vx

  x = - vdotnInv * (uYx * cpvx + uYy * cpvy + uYz * cpvz)
  y =   vdotnInv * (uXx * cpvx + uXy * cpvy + uXz * cpvz)
!
!
!     ...Shift to positive (x,y) pairs within the screen and rescale to [0-1]. Check, if
!        we want this proton recorded. Record also the detector number.
!
!
  x01 = (x + sHalf) * sInv
  y01 = (y + sHalf) * sInv

  hitScreen = (x01 >= 0.0 .and. x01 <= 1.0 .and. y01 >= 0.0 .and. y01 <= 1.0)
  recordProton = hitScreen .or. pi_recordOffscreenProtons
!
!
!     ...Record the proton (if at all) on the detector. Do so only if it goes through
!        the detector's associated pinhole (if present).
!
!
  if (recordProton) then

      pinholeRadius = pi_detectors (detector) % pinholeRadius

      if (pinholeRadius > 0.0) then

          hpx = pi_detectors (detector) % pinholeX - px
          hpy = pi_detectors (detector) % pinholeY - py
          hpz = pi_detectors (detector) % pinholeZ - pz

          hpvx = hpy * vz - hpz * vy
          hpvy = hpz * vx - hpx * vz
          hpvz = hpx * vy - hpy * vx

          hx = - vdotnInv * (uYx * hpvx + uYy * hpvy + uYz * hpvz)
          hy =   vdotnInv * (uXx * hpvx + uXy * hpvy + uXz * hpvz)

          recordProton = (hx*hx + hy*hy <= pinholeRadius * pinholeRadius)
      end if

      if (recordProton) then
          pi_screenProtonCount = pi_screenProtonCount + 1
          pi_screenProtons (SCREEN_POSX, pi_screenProtonCount) = x01
          pi_screenProtons (SCREEN_POSY, pi_screenProtonCount) = y01
          pi_screenProtons (SCREEN_DETC, pi_screenProtonCount) = real (detector)
          pi_screenProtons (SCREEN_DGJV, pi_screenProtonCount) = Jv
          pi_screenProtons (SCREEN_DGKX, pi_screenProtonCount) = Kx
          pi_screenProtons (SCREEN_DGKY, pi_screenProtonCount) = Ky
          pi_screenProtons (SCREEN_DGKZ, pi_screenProtonCount) = Kz
      end if

  end if
!
!
!     ...Optionally return info about screen hit and global screen coordinates.
!
!
  if (present (onScreen)) onScreen = hitScreen .and. recordProton
  if (present (sx)) sx = cx + x * uXx + y * uYx
  if (present (sy)) sy = cy + x * uXy + y * uYy
  if (present (sz)) sz = cz + x * uXz + y * uYz
!
!
!     ...Ready!
!
!
  return
end subroutine pi_recordProtonOnScreen
