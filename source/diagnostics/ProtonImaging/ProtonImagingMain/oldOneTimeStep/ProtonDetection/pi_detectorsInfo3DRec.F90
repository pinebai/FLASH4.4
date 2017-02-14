!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/oldOneTimeStep/ProtonDetection/pi_detectorsInfo3DRec
!!
!! NAME
!!
!!  pi_detectorsInfo3DRec
!!
!! SYNOPSIS
!!
!!  call pi_detectorsInfo3DRec ()
!!
!! DESCRIPTION
!!
!!  Generates information about the proton detectors placed in 3D rectangular (cartesian)
!!  space. In here all detector information is generated that can be generated at initialization.
!!  Currently it contains the following:
!!
!!         1) Calculate 1/2 and inverse of side length (for efficiency purposes).
!!
!!         2) Calculate the normal vector (if detector is aligned wrt to a proton beam).
!!
!!         3) Calculate the unit normal vector.
!!
!!         4) Calculate two orthogonal (x,y) axis unit vectors within the square screen
!!            plane of each detector, which will serve as a 2D cartesian local basis for
!!            the screen, located at the screen center.
!!
!!         5) Calculate the global coordinates of the four corners of each detector. The
!!            corners are labeled according to the direction they are located within the
!!            axis unit vectors when looking along the unit normal vector:
!!
!!                 right / left  label -> [+1,-1] direction of x-axis unit vectors
!!                 upper / lower label -> [+1,-1] direction of y-axis unit vectors
!!
!! ARGUMENTS
!!
!!  none
!!
!! NOTES
!!
!!***

subroutine pi_detectorsInfo3DRec ()

  use Driver_interface,    ONLY : Driver_abortFlash

  use Logfile_interface,   ONLY : Logfile_stamp

  use ProtonImaging_data,  ONLY : pi_badTiltingAxis,        &
                                  pi_beams,                 &
                                  pi_detectors,             &
                                  pi_normalizedTolerance,   &
                                  pi_numberOfDetectors,     &
                                  pi_orthogonalTolerance

  implicit none

#include "ProtonImaging.h"
#include "Flash.h"
#include "constants.h"

  character (len = 1) :: sideTiltingAxis

  integer :: alignWRTbeamNr
  integer :: detector

  real    :: centerX, centerY, centerZ
  real    :: capsuleX, capsuleY, capsuleZ
  real    :: cosPhi, sinPhi
  real    :: distance2beamCapsule
  real    :: normX, normY
  real    :: normalLengthInv
  real    :: nProj, nProjInv
  real    :: nx, ny, nz, nxny, nxnz, nynz
  real    :: orthoXY, orthoXn, orthoYn
  real    :: pd, px, py, pz
  real    :: sideLength, sideLengthHalf
  real    :: sideTiltingAngle
  real    :: targetX, targetY, targetZ
  real    :: uXx, uXy, uXz
  real    :: uYx, uYy, uYz
!
!
!     ...Loop over all detectors.
!
!
  do detector = 1, pi_numberOfDetectors

     centerX = pi_detectors (detector) % centerX
     centerY = pi_detectors (detector) % centerY
     centerZ = pi_detectors (detector) % centerZ
!
!
!     ...Calculate 1/2 and inverse of the detectors side length. This avoids having to
!        recompute it for each proton during screen detection.
!
!
     sideLength     = pi_detectors (detector) % sideLength
     sideLengthHalf = 0.5 * sideLength

     pi_detectors (detector) % sideLengthHalf = sideLengthHalf
     pi_detectors (detector) % sideLengthInv  = 1.0 / sideLength
!
!
!     ...If the detector should be aligned wrt a beam, calculate the unit normal vector
!        and the center position of the detector screen. Otherwise form the unit normal
!        vector from the normal vector.
!
!
     alignWRTbeamNr = pi_detectors (detector) % alignWRTbeamNr

     if (alignWRTbeamNr > 0) then

         capsuleX = pi_beams (alignWRTbeamNr) % capsuleX
         capsuleY = pi_beams (alignWRTbeamNr) % capsuleY
         capsuleZ = pi_beams (alignWRTbeamNr) % capsuleZ
         targetX  = pi_beams (alignWRTbeamNr) % targetX
         targetY  = pi_beams (alignWRTbeamNr) % targetY
         targetZ  = pi_beams (alignWRTbeamNr) % targetZ

         nx = targetX - capsuleX    ! detector normal vector pointing from beam capsule to target
         ny = targetY - capsuleY
         nz = targetZ - capsuleZ

         normalLengthInv = 1.0 / sqrt (nx * nx + ny * ny + nz * nz)

         nx = nx * normalLengthInv
         ny = ny * normalLengthInv
         nz = nz * normalLengthInv

         distance2beamCapsule = pi_detectors (detector) % distance2beamCapsule

         pi_detectors (detector) % centerX = capsuleX + nx * distance2beamCapsule
         pi_detectors (detector) % centerY = capsuleY + ny * distance2beamCapsule
         pi_detectors (detector) % centerZ = capsuleZ + nz * distance2beamCapsule

     else

         nx = pi_detectors (detector) % normalX
         ny = pi_detectors (detector) % normalY
         nz = pi_detectors (detector) % normalZ

         normalLengthInv = 1.0 / sqrt (nx * nx + ny * ny + nz * nz)

         nx = nx * normalLengthInv
         ny = ny * normalLengthInv
         nz = nz * normalLengthInv

     end if

     pi_detectors (detector) % normalX = nx
     pi_detectors (detector) % normalY = ny
     pi_detectors (detector) % normalZ = nz
!
!
!     ...If a pinhole is present, calculate the pinhole center coordinates.
!
!
     if (pi_detectors (detector) % pinholeRadius > 0.0) then

         pd = pi_detectors (detector) % pinholeDist2Det

         px = centerX - pd * nx     !
         py = centerY - pd * ny     ! minus sign ensures pinhole placement opposite of normal vector
         pz = centerZ - pd * nz     !

         pi_detectors (detector) % pinholeX = px
         pi_detectors (detector) % pinholeY = py
         pi_detectors (detector) % pinholeZ = pz

     end if
!
!
!     ...Set up 2 orthogonal unit vectors uX and uY (the 2D detector screen basis vectors) in
!        terms of the local x,y,z-coordinate system, whose origin coincides with center location
!        of the detector screen:
!
!
!                                    |\
!                                    | \
!                                    |  \
!                                    |   \
!                                    | uY|\
!                                    |   | \                C = detector center
!                                    |   |  \
!                      detector      |   |   |
!                                    |   C---|---> n        n = normal unit vector
!                       screen       \    \  |
!                                     \    \ |
!                                      \  uX\|        uX,uY,n = form a right hand rule
!                                       \    |                  orthogonal vector system:
!                                        \   |                      n = uY x uX
!                                         \  |
!                                          \ |
!                                           \|
!
!
!        The Y unit vector uY will be aligned and tilted with respect to the tilting axis
!        (one of the global X,Y,Z axis). The tilting angle A is defined as follows: 1) initially
!        (tilting angle = 0) the three vectors uY, n and the tilting axis are coplanar, 2) the
!        tilting angle is the clockwise angle rotation of uY from the tilting axis along the
!        normal vector n:
!
!
!                                 tilting
!                                  axis
!                                    |
!                                    |     uY           This view is along the tail T of
!                                    |     /            the normal vector.
!                                    |    /
!                                    | A /
!                                    |  /
!                                    | /
!                                    |/
!                                    T (n)
!
!        In case the tilting axis is parallel to the normal vector, the tilting angle becomes
!        undefined and another tilting axis must be chosen. The code below checks this and
!        informs the user.
!
!
     sideTiltingAngle = pi_detectors (detector) % sideTiltingAngle
     sideTiltingAxis  = pi_detectors (detector) % sideTiltingAxis

     sinPhi = sin (sideTiltingAngle)          ! tilting angle is in radians here
     cosPhi = cos (sideTiltingAngle)

     select case (sideTiltingAxis)

     case ('x')

       nxny  = nx * ny
       nxnz  = nx * nz
       nProj = sqrt (ny * ny + nz * nz)

       if (nProj == 0.0) then
           call Driver_abortFlash ("pi_detectorsInfo3DRec: Impossible detector tilting x-axis!")
       else if (nProj < pi_badTiltingAxis) then
           call Logfile_stamp ('Bad detector tilting axis, but will proceed...','[pi_detectorsInfo3DRec]')
       end if

       nProjInv = 1.0 / nProj

       uXx =  - nProj * sinPhi
       uXy = nProjInv * (  nxny * sinPhi  +  nz * cosPhi)
       uXz = nProjInv * (  nxnz * sinPhi  -  ny * cosPhi)

       uYx =    nProj * cosPhi
       uYy = nProjInv * (- nxny * cosPhi  +  nz * sinPhi)
       uYz = nProjInv * (- nxnz * cosPhi  -  ny * sinPhi)

     case ('y')

       nxny  = nx * ny
       nynz  = ny * nz
       nProj = sqrt (nx * nx + nz * nz)

       if (nProj == 0.0) then
           call Driver_abortFlash ("pi_detectorsInfo3DRec: Impossible detector tilting y-axis!")
       else if (nProj < pi_badTiltingAxis) then
           call Logfile_stamp ('Bad detector tilting axis, but will proceed...','[pi_detectorsInfo3DRec]')
       end if

       nProjInv = 1.0 / nProj

       uXx = nProjInv * (  nxny * sinPhi  -  nz * cosPhi)
       uXy =  - nProj * sinPhi
       uXz = nProjInv * (  nynz * sinPhi  +  nx * cosPhi)

       uYx = nProjInv * (- nxny * cosPhi  -  nz * sinPhi)
       uYy =    nProj * cosPhi
       uYz = nProjInv * (- nynz * cosPhi  +  nx * sinPhi)

     case ('z')

       nxnz  = nx * nz
       nynz  = ny * nz
       nProj = sqrt (nx * nx + ny * ny)

       if (nProj == 0.0) then
           call Driver_abortFlash ("pi_detectorsInfo3DRec: Impossible detector tilting z-axis!")
       else if (nProj < pi_badTiltingAxis) then
           call Logfile_stamp ('Bad detector tilting axis, but will proceed...','[pi_detectorsInfo3DRec]')
       end if

       nProjInv = 1.0 / nProj

       uXx = nProjInv * (  nxnz * sinPhi  +  ny * cosPhi)
       uXy = nProjInv * (  nynz * sinPhi  -  nx * cosPhi)
       uXz =  - nProj * sinPhi

       uYx = nProjInv * (- nxnz * cosPhi  +  ny * sinPhi)
       uYy = nProjInv * (- nynz * cosPhi  -  nx * sinPhi)
       uYz =    nProj * cosPhi

     case default

       call Driver_abortFlash ("pi_detectorsInfo3DRec: Invalid detector tilting axis!")

     end select
!
!
!     ...unit vector components are ready. Check for unit vector conditions and if ok
!        store into detectors array.
!
!
     normX   = uXx * uXx  +  uXy * uXy  +  uXz * uXz
     normY   = uYx * uYx  +  uYy * uYy  +  uYz * uYz
     orthoXY = uXx * uYx  +  uXy * uYy  +  uXz * uYz
     orthoXn = uXx *  nx  +  uXy *  ny  +  uXz *  nz
     orthoYn = uYx *  nx  +  uYy *  ny  +  uYz *  nz

     if (     (normX > 1.0 + pi_normalizedTolerance) &
         .or. (normX < 1.0 - pi_normalizedTolerance) ) then
          call Driver_abortFlash ("pi_detectorsInfo3DRec: Axis X unit vector not normalized!")
     end if

     if (     (normY > 1.0 + pi_normalizedTolerance) &
         .or. (normY < 1.0 - pi_normalizedTolerance) ) then
          call Driver_abortFlash ("pi_detectorsInfo3DRec: Axis Y unit vector not normalized!")
     end if

     if (     (orthoXY >   pi_orthogonalTolerance) &
         .or. (orthoXY < - pi_orthogonalTolerance) ) then
          call Driver_abortFlash ("pi_detectorsInfo3DRec: Axis X,Y unit vectors not orthogonal!")
     end if

     if (     (orthoXn >   pi_orthogonalTolerance) &
         .or. (orthoXn < - pi_orthogonalTolerance) ) then
          call Driver_abortFlash ("pi_detectorsInfo3DRec: Axis X + detector normal vectors not orthogonal!")
     end if

     if (     (orthoYn >   pi_orthogonalTolerance) &
         .or. (orthoYn < - pi_orthogonalTolerance) ) then
          call Driver_abortFlash ("pi_detectorsInfo3DRec: Axis Y + detector normal vectors not orthogonal!")
     end if

     pi_detectors (detector) % axisXunitX = uXx
     pi_detectors (detector) % axisXunitY = uXy
     pi_detectors (detector) % axisXunitZ = uXz
     pi_detectors (detector) % axisYunitX = uYx
     pi_detectors (detector) % axisYunitY = uYy
     pi_detectors (detector) % axisYunitZ = uYz
!
!
!     ...Calculate the global coordinates of the four corners. Using the center C global
!        X,Y,Z coordinates and the local screen x,y,z coordinates, each corner K has the
!        following global coordinate:
!
!                   K (X,Y,Z) = C (X,Y,Z) + S/2 * [+/- uX (x,y,z) +/- uY (x,y,z)]
!
!        where S is the screen's side length and uX,uY are the unit axis vectors defining
!        the screen area:
!
!
!                                   K4\
!                                    | \
!                                    |  \
!                                    |   \
!                                    | uY|\
!                                  S |   | \             C = detector center
!                                    |   |  \            S = screen side length
!                      detector      |   |   K1          n = normal unit vector
!                                    |   C---|---> n     K1 = upper right corner
!                       screen      K3    \  |           K2 = lower right corner
!                                     \    \ |           K3 = lower  left corner
!                                      \  uX\|           K4 = upper  left corner
!                                       \    |           uX,uY = screen unit vectors
!                                        \   |
!                                         \  |
!                                          \ |
!                                           \K2
!
!
     pi_detectors (detector) % cornerUpperRightX = centerX + sideLengthHalf * ( + uXx + uYx)
     pi_detectors (detector) % cornerUpperRightY = centerY + sideLengthHalf * ( + uXy + uYy)
     pi_detectors (detector) % cornerUpperRightZ = centerZ + sideLengthHalf * ( + uXz + uYz)
     pi_detectors (detector) % cornerLowerRightX = centerX + sideLengthHalf * ( + uXx - uYx)
     pi_detectors (detector) % cornerLowerRightY = centerY + sideLengthHalf * ( + uXy - uYy)
     pi_detectors (detector) % cornerLowerRightZ = centerZ + sideLengthHalf * ( + uXz - uYz)
     pi_detectors (detector) % cornerUpperLeftX  = centerX + sideLengthHalf * ( - uXx + uYx)
     pi_detectors (detector) % cornerUpperLeftY  = centerY + sideLengthHalf * ( - uXy + uYy)
     pi_detectors (detector) % cornerUpperLeftZ  = centerZ + sideLengthHalf * ( - uXz + uYz)
     pi_detectors (detector) % cornerLowerLeftX  = centerX + sideLengthHalf * ( - uXx - uYx)
     pi_detectors (detector) % cornerLowerLeftY  = centerY + sideLengthHalf * ( - uXy - uYy)
     pi_detectors (detector) % cornerLowerLeftZ  = centerZ + sideLengthHalf * ( - uXz - uYz)

  end do
!
!
!     ...Ready!
!
!
  return
end subroutine pi_detectorsInfo3DRec
