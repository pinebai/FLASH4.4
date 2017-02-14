!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/ProtonBeams/pi_beamsInfo3DRec
!!
!! NAME
!!
!!  pi_beamsInfo3DRec
!!
!! SYNOPSIS
!!
!!  call pi_beamsInfo3DRec ()
!!
!! DESCRIPTION
!!
!!  Generates information about the proton beams for those domain geometries consisting
!!  formally of 3D rectangular grids (cartesian). In here all beam information is generated
!!  that can be generated at initialization. Currently it contains the following:
!!
!!         1) Convert proton energy from MeV to initial proton speed in units
!!            of cm/s.
!!
!!         2) Calculate the capsule center to target center distance.
!!
!!         3) Calculate the expected number of protons per capsule grain and, if
!!            necessary, re-adjust total number of protons per beam.
!!
!!         4) Calculate three orthogonal unit vectors, such that one unit vector lays
!!            on the capsule/target center line with direction from capsule -> target.
!!            The remaining 2 unit vectors are then in a plane orthogonal to the
!!            capsule/target center line.
!!
!! ARGUMENTS
!!
!!  none
!!
!! NOTES
!!
!!***

subroutine pi_beamsInfo3DRec ()

  use Driver_interface,    ONLY : Driver_abortFlash

  use Logfile_interface,   ONLY : Logfile_stamp

  use pi_interface,        ONLY : pi_capsuleGrainIndices2xyz, &
                                  pi_capsuleNextGrainIndices, &
                                  pi_capsuleTotalGrainCount

  use ProtonImaging_data,  ONLY : pi_beams,                 &
                                  pi_MeV2erg,               &
                                  pi_normalizedTolerance,   &
                                  pi_numberOfBeams,         &
                                  pi_orthogonalTolerance,   &
                                  pi_protonMass,            &
                                  pi_speedOfLight,          &
                                  pi_speedOfLightSquared,   &
                                  pi_totalProtons2Launch

  implicit none

#include "ProtonImaging.h"
#include "Flash.h"
#include "constants.h"

  logical :: planeXY, planeXZ
  logical :: valid

  integer :: beam
  integer :: capsuleGrainLevel
  integer :: capsuleNumberOfGrains
  integer :: i,j,k
  integer :: numberOfProtons
  integer :: numberOfProtonsPerGrain

  real    :: apertureAngle
  real    :: beamLength, beamLengthInv
  real    :: bProj, bProjInv
  real    :: bx, by, bz
  real    :: capsuleRadius
  real    :: capsuleX, capsuleY, capsuleZ
  real    :: E, mc2, v
  real    :: grainSize
  real    :: norm1, norm2, norm3
  real    :: ortho12, ortho13, ortho23
  real    :: protonEnergy
  real    :: targetRadius
  real    :: targetX, targetY, targetZ
  real    :: u1x, u1y, u1z
  real    :: u2x, u2y, u2z
  real    :: u3x, u3y, u3z
  real    :: x, y, z
  real    :: zerox, zeroy, zeroz
!
!
!     ...Loop over all beams.
!
!
  do beam = 1, pi_numberOfBeams

     apertureAngle       = pi_beams (beam) % apertureAngle
     capsuleGrainLevel   = pi_beams (beam) % capsuleGrainLevel
     capsuleRadius       = pi_beams (beam) % capsuleRadius
     capsuleX            = pi_beams (beam) % capsuleX
     capsuleY            = pi_beams (beam) % capsuleY
     capsuleZ            = pi_beams (beam) % capsuleZ
     numberOfProtons     = pi_beams (beam) % numberOfProtons
     protonEnergy        = pi_beams (beam) % protonEnergy
     targetRadius        = pi_beams (beam) % targetRadius
     targetX             = pi_beams (beam) % targetX
     targetY             = pi_beams (beam) % targetY
     targetZ             = pi_beams (beam) % targetZ
!
!
!     ...Set the dimensionality indicator of the beam.
!
!
     pi_beams (beam) % dimensionality = 3
!
!
!     ...Calculate the initital proton speed (in units of cm/s) from the supplied
!        proton energy in MeV.
!
!
     E   = protonEnergy * pi_MeV2erg
     mc2 = pi_protonMass * pi_speedOfLightSquared
     v   = sqrt (1.0 - (mc2 / (E + mc2)) ** 2)     ! in units of light speed

     pi_beams (beam) % initialProtonSpeed = v * pi_speedOfLight
!
!
!     ...Calculate distance from capsule center to target and store in beams array.
!
!
     bx = targetX - capsuleX        ! beam vector pointing from capsule to target
     by = targetY - capsuleY
     bz = targetZ - capsuleZ

     beamLength    = sqrt (bx * bx + by * by + bz * bz)
     beamLengthInv = 1.0 / beamLength

     pi_beams (beam) % distanceCapsule2Target = beamLength
!
!
!     ...Set the random number generator seed for the beam. Currently this
!        is set equal to the beam number.
!
!
     pi_beams (beam) % randomNumberSeed = beam
!
!
!     ...If not set as a runtime parameter, calculate and store the target radius length
!        using the supplied aperture angle.
!
!
     if (targetRadius < 0.0) then
         pi_beams (beam) % targetRadius = tan (apertureAngle * 0.5) * beamLength
     end if
!
!
!     ...Set up 3 orthogonal unit vectors with the 3rd unit vector u3 along the beam
!        center line and pointing from capsule to target:
!
!
!                                    T
!                                   /
!                                  /
!                                 /
!                                /
!                          u1   /            u1 x u2 = u3
!                           |  u3
!                           | /
!                           |/
!                           C------u2
!
!
!        The 1st unit vector u1 will be aligned with one of the global X,Y,Z axis, the
!        choice of the axis depending on the magnitude of the beam vector projection onto
!        the XY,XZ,YZ planes:
!
!               XY-plane projection of beam vector largest -> align along Z axis
!               XZ-plane projection of beam vector largest -> align along Y axis
!               YZ-plane projection of beam vector largest -> align along X axis
!
!               align along Z axis -> means: u1 torsion angle with Z axis is zero
!               align along Y axis -> means: u1 torsion angle with Y axis is zero
!               align along X axis -> means: u1 torsion angle with X axis is zero
!
!        The 2nd unit vector u2 will be placed at right angles to u1 in clockwise rotation
!        when looking along the beam line from capsule -> target. The 3D unit system u1,u2,u3
!        will form a right handed coordinate system in which u1 x u2 = u3.
!
!
     u3x = bx * beamLengthInv
     u3y = by * beamLengthInv
     u3z = bz * beamLengthInv

     planeXY = (abs (bx) >= abs (by)) .and. (abs (by) >= abs (bz))      ! |bz| is smallest
     planeXZ = (abs (bx) >= abs (by)) .and. (abs (by) <  abs (bz))      ! |by| is smallest

     if (planeXY) then

         bProj    = sqrt (bx * bx + by * by)
         bProjInv = 1.0 / bProj

         u1x = - bx * bz * bProjInv * beamLengthInv
         u1y = - by * bz * bProjInv * beamLengthInv
         u1z =   bProj * beamLengthInv

         u2x =   by * bProjInv
         u2y = - bx * bProjInv
         u2z = 0.0

     else if (planeXZ) then

         bProj    = sqrt (bx * bx + bz * bz)
         bProjInv = 1.0 / bProj

         u1x = - bx * by * bProjInv * beamLengthInv
         u1y =   bProj * beamLengthInv
         u1z = - by * bz * bProjInv * beamLengthInv

         u2x = - bz * bProjInv
         u2y = 0.0
         u2z =   bx * bProjInv

     else

         bProj    = sqrt (by * by + bz * bz)
         bProjInv = 1.0 / bProj

         u1x =   bProj * beamLengthInv
         u1y = - bx * by * bProjInv * beamLengthInv
         u1z = - bx * bz * bProjInv * beamLengthInv

         u2x = 0.0
         u2y =   bz * bProjInv
         u2z = - by * bProjInv

     end if
!
!
!     ...unit vector components are ready. Check for unit vector conditions and if ok
!        store into beams array.
!
!
     norm1   = u1x * u1x  +  u1y * u1y  +  u1z * u1z
     norm2   = u2x * u2x  +  u2y * u2y  +  u2z * u2z
     norm3   = u3x * u3x  +  u3y * u3y  +  u3z * u3z
     ortho12 = u1x * u2x  +  u1y * u2y  +  u1z * u2z
     ortho13 = u1x * u3x  +  u1y * u3y  +  u1z * u3z
     ortho23 = u2x * u3x  +  u2y * u3y  +  u2z * u3z
     zerox   = u3x - (u1y * u2z  -  u1z * u2y)
     zeroy   = u3y - (u1z * u2x  -  u1x * u2z)
     zeroz   = u3z - (u1x * u2y  -  u1y * u2x)

     if (     (norm1 > 1.0 + pi_normalizedTolerance) &
         .or. (norm1 < 1.0 - pi_normalizedTolerance) ) then
          call Driver_abortFlash ("pi_beamsInfo3DRec: 1st unit vector not normalized!")
     end if

     if (     (norm2 > 1.0 + pi_normalizedTolerance) &
         .or. (norm2 < 1.0 - pi_normalizedTolerance) ) then
          call Driver_abortFlash ("pi_beamsInfo3DRec: 2nd unit vector not normalized!")
     end if

     if (     (norm3 > 1.0 + pi_normalizedTolerance) &
         .or. (norm3 < 1.0 - pi_normalizedTolerance) ) then
          call Driver_abortFlash ("pi_beamsInfo3DRec: 3rd unit vector not normalized!")
     end if

     if (     (ortho12 >   pi_orthogonalTolerance) &
         .or. (ortho12 < - pi_orthogonalTolerance) ) then
          call Driver_abortFlash ("pi_beamsInfo3DRec: 1st + 2nd unit vectors not orthogonal!")
     end if

     if (     (ortho13 >   pi_orthogonalTolerance) &
         .or. (ortho13 < - pi_orthogonalTolerance) ) then
          call Driver_abortFlash ("pi_beamsInfo3DRec: 1st + 3rd unit vectors not orthogonal!")
     end if

     if (     (ortho23 >   pi_orthogonalTolerance) &
         .or. (ortho23 < - pi_orthogonalTolerance) ) then
          call Driver_abortFlash ("pi_beamsInfo3DRec: 2nd + 3rd unit vectors not orthogonal!")
     end if

     if (     abs (zerox) > pi_orthogonalTolerance &
         .or. abs (zeroy) > pi_orthogonalTolerance &
         .or. abs (zeroz) > pi_orthogonalTolerance ) then
          call Driver_abortFlash ("pi_beamsInfo3DRec: u1 x u2 is not equal to u3!")
     end if

     pi_beams (beam) % axisUnit1X = u1x
     pi_beams (beam) % axisUnit1Y = u1y
     pi_beams (beam) % axisUnit1Z = u1z
     pi_beams (beam) % axisUnit2X = u2x
     pi_beams (beam) % axisUnit2Y = u2y
     pi_beams (beam) % axisUnit2Z = u2z
     pi_beams (beam) % axisUnit3X = u3x
     pi_beams (beam) % axisUnit3Y = u3y
     pi_beams (beam) % axisUnit3Z = u3z
!
!
!     ...Calculate the total number of grains in capsule, the expected number of protons
!        per capsule grain and, if necessary, re-adjust total number of protons per beam.
!        Find and set the first valid capsule grain indices and prepare for launching
!        the first set of protons from that grain.
!
!
     if (capsuleGrainLevel == 0) then

         capsuleNumberOfGrains   = numberOfProtons
         numberOfProtonsPerGrain = 1

     else if (capsuleGrainLevel == 1) then

         i = 0
         j = 0
         k = 0

         capsuleNumberOfGrains   = 1
         numberOfProtonsPerGrain = numberOfProtons

         pi_beams (beam) % capsuleGrainIndexI = i
         pi_beams (beam) % capsuleGrainIndexJ = j
         pi_beams (beam) % capsuleGrainIndexK = k
         pi_beams (beam) % capsuleGrainLocalX = 0.0
         pi_beams (beam) % capsuleGrainLocalY = 0.0
         pi_beams (beam) % capsuleGrainLocalZ = 0.0
         pi_beams (beam) % capsuleGrainSize   = capsuleRadius

     else
         capsuleNumberOfGrains = pi_capsuleTotalGrainCount (capsuleGrainLevel)

         if (capsuleNumberOfGrains < 0) then
             call Driver_abortFlash ("pi_beamsInfo3DRec: Capsule number of grains overflow!")
         else if (capsuleNumberOfGrains == 0) then
             call Driver_abortFlash ("pi_beamsInfo3DRec: Capsule number of grains = 0!")
         end if

         numberOfProtonsPerGrain = ceiling (real (numberOfProtons) / real (capsuleNumberOfGrains))
         numberOfProtons         = numberOfProtonsPerGrain * capsuleNumberOfGrains

         i = - capsuleGrainLevel
         j = - capsuleGrainLevel
         k = - capsuleGrainLevel

         call pi_capsuleNextGrainIndices (capsuleGrainLevel, i,j,k,valid)

         if (.not.valid) then
             call Driver_abortFlash ("pi_beamsInfo3DRec: No initial capsule grain indices!")
         end if

         grainSize = capsuleRadius / real (capsuleGrainLevel)

         call pi_capsuleGrainIndices2xyz (i,j,k,grainSize,   x,y,z)   ! capsule grain local x,y,z

         pi_beams (beam) % capsuleGrainIndexI = i
         pi_beams (beam) % capsuleGrainIndexJ = j
         pi_beams (beam) % capsuleGrainIndexK = k
         pi_beams (beam) % capsuleGrainLocalX = x
         pi_beams (beam) % capsuleGrainLocalY = y
         pi_beams (beam) % capsuleGrainLocalZ = z
         pi_beams (beam) % capsuleGrainSize   = grainSize

     end if

     pi_beams (beam) % capsuleGrainProtons     = numberOfProtonsPerGrain
     pi_beams (beam) % capsuleNumberOfGrains   = capsuleNumberOfGrains
     pi_beams (beam) % numberOfProtons         = numberOfProtons
     pi_beams (beam) % numberOfProtonsPerGrain = numberOfProtonsPerGrain

  end do
!
!
!     ...Calculate the total number of protons that will be launched from all beams.
!
!
  pi_totalProtons2Launch = sum (pi_beams (1:pi_numberOfBeams) % numberOfProtons)
!
!
!     ...Ready!
!
!
  return
end subroutine pi_beamsInfo3DRec
