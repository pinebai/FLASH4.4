!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserBeams/ed_beamsInfo3DRec
!!
!! NAME
!!
!!  ed_beamsInfo3DRec
!!
!! SYNOPSIS
!!
!!  call ed_beamsInfo3DRec ()
!!
!! DESCRIPTION
!!
!!  Generates information about the beams for those domain geometries consisting formally of
!!  3D rectangular grids (cartesian). In here all beam information is generated that can be
!!  generated at initialization. Currently it contains the following:
!!
!!         1) Calculate the lens to target distance and incorporate some
!!            additional info (time duration) regarding the pulses associated
!!            with each beam.
!!
!!         2) Calculate basic information about the planar beam ray grid.
!!
!!         3) Calculate the unit vectors associated with both the semiaxes
!!            in the local target coordinate system (i.e. with origin on
!!            the center of the ellipse/rectangle).
!!
!! ARGUMENTS
!!
!! NOTES
!!
!!***

subroutine ed_beamsInfo3DRec ()

  use Driver_interface,       ONLY : Driver_abortFlash

  use Logfile_interface,      ONLY : Logfile_stamp

  use ed_interface,           ONLY : ed_beam3DGridSetupDelta,       &
                                     ed_beam3DGridSetupRadial,      &
                                     ed_beam3DGridSetupRecBeam,     &
                                     ed_beam3DGridSetupSquare,      &
                                     ed_beam3DGridSetupStatistical, &
                                     ed_beam3DGridWeightDelta,      &
                                     ed_beam3DGridWeightRadial,     &
                                     ed_beam3DGridWeightRecBeam,    &
                                     ed_beam3DGridWeightSquare

  use EnergyDeposition_data,  ONLY : ed_badTorsionAxis,        &
                                     ed_beams,                 &
                                     ed_normalizedTolerance,   &
                                     ed_numberOfBeams,         &
                                     ed_orthogonalTolerance,   &
                                     ed_pulseNumberOfSections, &
                                     ed_pulses

  implicit none

#include "EnergyDeposition.h"
#include "Flash.h"
#include "constants.h"

  character (len = BEAM_STRING_LENGTH) :: functionType
  character (len = BEAM_STRING_LENGTH) :: gridType
  character (len = 1                 ) :: torsionAxis

  logical :: seedIncrement, seedInitialize

  integer :: beam
  integer :: lastSection
  integer :: nRays
  integer :: nTics1stDim, nTics2ndDim
  integer :: pulseNumber
  integer :: seed, seedMaximum, seedStepping

  real    :: beamLength
  real    :: beginTime, endTime
  real    :: bProj, bProjInv
  real    :: bsizeInv
  real    :: bx, by, bz
  real    :: bxby, bxbz, bybz
  real    :: delta1stDim, delta2ndDim
  real    :: firstTic1stDim, firstTic2ndDim
  real    :: gaussianCenterMajor, gaussianCenterMinor
  real    :: gaussianExponent
  real    :: gaussianRadiusMajor, gaussianRadiusMinor
  real    :: gridWeight
  real    :: lensSemiAxisMajor, lensSemiAxisMinor
  real    :: lensX, lensY, lensZ
  real    :: norm1, norm2
  real    :: ortho12, ortho1b, ortho2b
  real    :: sinPhi, cosPhi
  real    :: target2LensMagnify
  real    :: targetSemiAxisMajor, targetSemiAxisMinor
  real    :: targetX, targetY, targetZ
  real    :: torsionAngle
  real    :: u1x, u1y, u1z
  real    :: u2x, u2y, u2z
  real    :: x,y,z
!
!
!     ...Loop over all beams.
!
!
  do beam = 1, ed_numberOfBeams

     delta1stDim         = ed_beams (beam) % gridDelta1stDim
     delta2ndDim         = ed_beams (beam) % gridDelta2ndDim
     functionType        = ed_beams (beam) % crossSectionFunctionType
     gaussianCenterMajor = ed_beams (beam) % gaussianCenterMajor
     gaussianCenterMinor = ed_beams (beam) % gaussianCenterMinor
     gaussianExponent    = ed_beams (beam) % gaussianExponent
     gaussianRadiusMajor = ed_beams (beam) % gaussianRadiusMajor
     gaussianRadiusMinor = ed_beams (beam) % gaussianRadiusMinor
     gridType            = ed_beams (beam) % gridType
     lensSemiAxisMajor   = ed_beams (beam) % lensSemiAxisMajor
     lensX               = ed_beams (beam) % lensX
     lensY               = ed_beams (beam) % lensY
     lensZ               = ed_beams (beam) % lensZ
     nRays               = ed_beams (beam) % numberOfRays
     nTics1stDim         = ed_beams (beam) % gridnTics1stDim
     nTics2ndDim         = ed_beams (beam) % gridnTics2ndDim
     pulseNumber         = ed_beams (beam) % pulseNumber
     targetSemiAxisMajor = ed_beams (beam) % targetSemiAxisMajor
     targetSemiAxisMinor = ed_beams (beam) % targetSemiAxisMinor
     targetX             = ed_beams (beam) % targetX
     targetY             = ed_beams (beam) % targetY
     targetZ             = ed_beams (beam) % targetZ
     torsionAngle        = ed_beams (beam) % semiAxisMajorTorsionAngle
     torsionAxis         = ed_beams (beam) % semiAxisMajorTorsionAxis
!
!
!     ...Set the dimensionality indicator of the beam.
!
!
     ed_beams (beam) % dimensionality = 3
!
!
!     ...Calculate distance from lens to target and pulse time bounds for each beam.
!
!
     x = targetX - lensX
     y = targetY - lensY
     z = targetZ - lensZ

     beamLength = sqrt (x * x + y * y + z * z)

     lastSection = ed_pulseNumberOfSections (pulseNumber)
     beginTime   = ed_pulses (1,             pulseNumber) % pulseTime
     endTime     = ed_pulses (lastSection,   pulseNumber) % pulseTime

     ed_beams (beam) % distanceLens2Target = beamLength
     ed_beams (beam) % pulseStartingTime   = beginTime
     ed_beams (beam) % pulseEndingTime     = endTime
!
!
!     ...Set up the beam's grid and store the info.
!
!
     if (gridType == 'delta2D') then

         call ed_beam3DGridSetupDelta   (targetSemiAxisMajor,              &
                                         targetSemiAxisMinor,              &
                                         delta1stDim,                      &
                                         delta2ndDim,                      &
                                                           nRays,          &
                                                           nTics1stDim,    &
                                                           nTics2ndDim,    &
                                                           firstTic1stDim, &
                                                           firstTic2ndDim  )

         call ed_beam3DGridWeightDelta  (targetSemiAxisMajor,              &
                                         targetSemiAxisMinor,              &
                                         functionType,                     &
                                         gaussianExponent,                 &
                                         gaussianRadiusMajor,              &
                                         gaussianRadiusMinor,              &
                                         gaussianCenterMajor,              &
                                         gaussianCenterMinor,              &
                                         nTics1stDim,                      &
                                         nTics2ndDim,                      &
                                         delta1stDim,                      &
                                         delta2ndDim,                      &
                                         firstTic1stDim,                   &
                                         firstTic2ndDim,                   &
                                         nRays,                            &
                                                                gridWeight )

         ed_beams (beam) % numberOfRays       = nRays
         ed_beams (beam) % gridnTics1stDim    = nTics1stDim
         ed_beams (beam) % gridnTics2ndDim    = nTics2ndDim
         ed_beams (beam) % gridFirstTic1stDim = firstTic1stDim
         ed_beams (beam) % gridFirstTic2ndDim = firstTic2ndDim
         ed_beams (beam) % gridWeight         = gridWeight

     else if (gridType == 'square2D') then

         call ed_beam3DGridSetupSquare  (targetSemiAxisMajor,              &
                                         targetSemiAxisMinor,              &
                                                              nRays,       &
                                                              nTics1stDim, &
                                                              nTics2ndDim, &
                                                              delta1stDim  )

         call ed_beam3DGridWeightSquare (targetSemiAxisMajor,              &
                                         targetSemiAxisMinor,              &
                                         functionType,                     &
                                         gaussianExponent,                 &
                                         gaussianRadiusMajor,              &
                                         gaussianRadiusMinor,              &
                                         gaussianCenterMajor,              &
                                         gaussianCenterMinor,              &
                                         nTics1stDim,                      &
                                         nTics2ndDim,                      &
                                         delta1stDim,                      &
                                         nRays,                            &
                                                                gridWeight )

         ed_beams (beam) % numberOfRays       = nRays
         ed_beams (beam) % gridnTics1stDim    = nTics1stDim
         ed_beams (beam) % gridnTics2ndDim    = nTics2ndDim
         ed_beams (beam) % gridDelta1stDim    = delta1stDim
         ed_beams (beam) % gridDelta2ndDim    = delta1stDim      ! same delta for 2nd dimension
         ed_beams (beam) % gridFirstTic1stDim = 0.0              ! no first tic displacement
         ed_beams (beam) % gridFirstTic2ndDim = 0.0              ! no first tic displacement
         ed_beams (beam) % gridWeight         = gridWeight

     else if (gridType == 'radial2D') then

         call ed_beam3DGridSetupRadial  (                  nRays,          &
                                                           nTics1stDim,    &
                                                           nTics2ndDim,    &
                                                           delta1stDim,    &
                                                           delta2ndDim,    &
                                                           firstTic1stDim, &
                                                           firstTic2ndDim  )

         call ed_beam3DGridWeightRadial (targetSemiAxisMajor,              &
                                         targetSemiAxisMinor,              &
                                         functionType,                     &
                                         gaussianExponent,                 &
                                         gaussianRadiusMajor,              &
                                         gaussianRadiusMinor,              &
                                         gaussianCenterMajor,              &
                                         gaussianCenterMinor,              &
                                         nTics1stDim,                      &
                                         nTics2ndDim,                      &
                                         delta1stDim,                      &
                                         delta2ndDim,                      &
                                         firstTic1stDim,                   &
                                         firstTic2ndDim,                   &
                                         nRays,                            &
                                                                gridWeight )

         ed_beams (beam) % numberOfRays       = nRays
         ed_beams (beam) % gridnTics1stDim    = nTics1stDim
         ed_beams (beam) % gridnTics2ndDim    = nTics2ndDim
         ed_beams (beam) % gridDelta1stDim    = delta1stDim
         ed_beams (beam) % gridDelta2ndDim    = delta2ndDim
         ed_beams (beam) % gridFirstTic1stDim = firstTic1stDim
         ed_beams (beam) % gridFirstTic2ndDim = firstTic2ndDim
         ed_beams (beam) % gridWeight         = gridWeight

     else if (gridType == 'rectangular2D') then

         call ed_beam3DGridSetupRecBeam  (targetSemiAxisMajor,              &
                                          targetSemiAxisMinor,              &
                                                               nRays,       &
                                                               nTics1stDim, &
                                                               nTics2ndDim, &
                                                               delta1stDim, &
                                                               delta2ndDim  )

         call ed_beam3DGridWeightRecBeam (functionType,                     &
                                          gaussianExponent,                 &
                                          gaussianRadiusMajor,              &
                                          gaussianRadiusMinor,              &
                                          gaussianCenterMajor,              &
                                          gaussianCenterMinor,              &
                                          nTics1stDim,                      &
                                          nTics2ndDim,                      &
                                          delta1stDim,                      &
                                          delta2ndDim,                      &
                                          nRays,                            &
                                                                 gridWeight )

         ed_beams (beam) % numberOfRays       = nRays
         ed_beams (beam) % gridnTics1stDim    = nTics1stDim
         ed_beams (beam) % gridnTics2ndDim    = nTics2ndDim
         ed_beams (beam) % gridDelta1stDim    = delta1stDim
         ed_beams (beam) % gridDelta2ndDim    = delta2ndDim
         ed_beams (beam) % gridFirstTic1stDim = 0.0              ! no first tic displacement
         ed_beams (beam) % gridFirstTic2ndDim = 0.0              ! no first tic displacement
         ed_beams (beam) % gridWeight         = gridWeight

     else if (gridType == 'statistical2D') then

         seedInitialize = .true.
         seedIncrement  = .false.

         call ed_beam3DGridSetupStatistical (seedInitialize,                &
                                             seedIncrement,                 &
                                                              seedMaximum,  &
                                                              seedStepping, &
                                                              seed          )

         ed_beams (beam) % gridSeed           = seed
         ed_beams (beam) % gridSeedMaximum    = seedMaximum
         ed_beams (beam) % gridSeedStepping   = seedStepping

     else
         call Driver_abortFlash ('[ed_beamsInfo3DRec] ERROR: unknown 3D beam grid type')
     end if
!
!
!     ...Calculate the magnification factor for target -> lens semiaxes. This magnification
!        factor will be needed to establish the beam grid data at the lens, since all the
!        beam grid data that has been stored is based on the target. Calculate and store
!        also the missing lens minor semiaxis length.
!
!
     target2LensMagnify = lensSemiAxisMajor  / targetSemiAxisMajor
     lensSemiAxisMinor  = target2LensMagnify * targetSemiAxisMinor

     ed_beams (beam) % lensSemiAxisMinor        = lensSemiAxisMinor
     ed_beams (beam) % target2LensMagnification = target2LensMagnify
!
!
!     ...Set up unit vectors in terms of the local target x,y,z-coordinate system (origin
!        coincides with origin of the two semiaxes inside elliptical target area).
!
!
     bx = lensX - targetX         ! beam vector pointing from target to lens
     by = lensY - targetY
     bz = lensZ - targetZ

     bSizeInv = 1.0 / beamLength

     sinPhi = sin (torsionAngle)          ! torsion angle is in radians here
     cosPhi = cos (torsionAngle)

     select case (torsionAxis)

     case ('x')

       bxby  = bx * by
       bxbz  = bx * bz
       bProj = sqrt (by * by + bz * bz)

       if (bProj == 0.0) then
           call Driver_abortFlash ("ed_beamsInfo3DRec: Impossible beam torsion z-axis!")
       end if

       if (bProj * bSizeInv < ed_badTorsionAxis) then
           call Logfile_stamp ('Bad beam torsion axis, but will proceed...','[ed_beamsInfo3DRec]')
       end if

       bProjInv = 1.0 / bProj

       u1x =    bProj * bSizeInv * cosPhi
       u1y = bProjInv * (- bxby * bSizeInv * cosPhi  -  bz * sinPhi)
       u1z = bProjInv * (- bxbz * bSizeInv * cosPhi  +  by * sinPhi)

       u2x =  - bProj * bSizeInv * sinPhi
       u2y = bProjInv * (  bxby * bSizeInv * sinPhi  -  bz * cosPhi)
       u2z = bProjInv * (  bxbz * bSizeInv * sinPhi  +  by * cosPhi)

     case ('y')

       bxby  = bx * by
       bybz  = by * bz
       bProj = sqrt (bx * bx + bz * bz)

       if (bProj == 0.0) then
           call Driver_abortFlash ("ed_beamsInfo3DRec: Impossible beam torsion y-axis!")
       else if (bProj * bSizeInv < ed_badTorsionAxis) then
           call Logfile_stamp ('Bad beam torsion axis, but will proceed...','[ed_beamsInfo3DRec]')
       end if

       bProjInv = 1.0 / bProj

       u1x = bProjInv * (- bxby * bSizeInv * cosPhi  +  bz * sinPhi)
       u1y =    bProj * bSizeInv * cosPhi
       u1z = bProjInv * (- bybz * bSizeInv * cosPhi  -  bx * sinPhi)

       u2x = bProjInv * (  bxby * bSizeInv * sinPhi  +  bz * cosPhi)
       u2y =  - bProj * bSizeInv * sinPhi
       u2z = bProjInv * (  bybz * bSizeInv * sinPhi  -  bx * cosPhi)

     case ('z')

       bxbz  = bx * bz
       bybz  = by * bz
       bProj = sqrt (bx * bx + by * by)

       if (bProj == 0.0) then
           call Driver_abortFlash ("ed_beamsInfo3DRec: Impossible beam torsion z-axis!")
       else if (bProj * bSizeInv < ed_badTorsionAxis) then
           call Logfile_stamp ('Bad beam torsion axis, but will proceed...','[ed_beamsInfo3DRec]')
       end if

       bProjInv = 1.0 / bProj

       u1x = bProjInv * (- bxbz * bSizeInv * cosPhi  -  by * sinPhi)
       u1y = bProjInv * (- bybz * bSizeInv * cosPhi  +  bx * sinPhi)
       u1z =    bProj * bSizeInv * cosPhi

       u2x = bProjInv * (  bxbz * bSizeInv * sinPhi  -  by * cosPhi)
       u2y = bProjInv * (  bybz * bSizeInv * sinPhi  +  bx * cosPhi)
       u2z =  - bProj * bSizeInv * sinPhi

     case default

       call Driver_abortFlash ("ed_beamsInfo3DRec: Invalid beam torsion axis!")

     end select
!
!
!     ...unit vector components are ready. Check for unit vector conditions and if ok
!        store into beams array.
!
!
     norm1   = u1x * u1x  +  u1y * u1y  +  u1z * u1z
     norm2   = u2x * u2x  +  u2y * u2y  +  u2z * u2z
     ortho12 = u1x * u2x  +  u1y * u2y  +  u1z * u2z
     ortho1b = u1x *  bx  +  u1y *  by  +  u1z *  bz
     ortho2b = u2x *  bx  +  u2y *  by  +  u2z *  bz

     if (     (norm1 > 1.0 + ed_normalizedTolerance) &
         .or. (norm1 < 1.0 - ed_normalizedTolerance) ) then
          call Driver_abortFlash ("ed_beamsInfo3DRec: 1st elliptical unit vector not normalized!")
     end if

     if (     (norm2 > 1.0 + ed_normalizedTolerance) &
         .or. (norm2 < 1.0 - ed_normalizedTolerance) ) then
          call Driver_abortFlash ("ed_beamsInfo3DRec: 2nd elliptical unit vector not normalized!")
     end if

     if (     (ortho12 >   ed_orthogonalTolerance) &
         .or. (ortho12 < - ed_orthogonalTolerance) ) then
          call Driver_abortFlash ("ed_beamsInfo3DRec: 1st + 2nd elliptical unit vectors not orthogonal!")
     end if

     if (     (ortho1b >   ed_orthogonalTolerance) &
         .or. (ortho1b < - ed_orthogonalTolerance) ) then
          call Driver_abortFlash ("ed_beamsInfo3DRec: 1st elliptical unit + beam vectors not orthogonal!")
     end if

     if (     (ortho2b >   ed_orthogonalTolerance) &
         .or. (ortho2b < - ed_orthogonalTolerance) ) then
          call Driver_abortFlash ("ed_beamsInfo3DRec: 2nd elliptical unit + beam vectors not orthogonal!")
     end if

     ed_beams (beam) % semiAxisUnitMajorX = u1x
     ed_beams (beam) % semiAxisUnitMajorY = u1y
     ed_beams (beam) % semiAxisUnitMajorZ = u1z
     ed_beams (beam) % semiAxisUnitMinorX = u2x
     ed_beams (beam) % semiAxisUnitMinorY = u2y
     ed_beams (beam) % semiAxisUnitMinorZ = u2z

  end do
!
!
!     ...Ready!
!
!
  return
end subroutine ed_beamsInfo3DRec
