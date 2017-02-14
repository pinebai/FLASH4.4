!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserBeams/ed_beamsInfo2DRec
!!
!! NAME
!!
!!  ed_beamsInfo2DRec
!!
!! SYNOPSIS
!!
!!  call ed_beamsInfo2DRec ()
!!
!! DESCRIPTION
!!
!!  Generates information about the beams for those geometries consisting formally of 2D
!!  rectangular grids (cartesian + cylindrical). In here all beam information is generated
!!  that can be generated at initialization. Currently it contains the following:
!!
!!         1) Calculate the lens to target distance and incorporate some
!!            additional info (time duration) regarding the pulses associated
!!            with each beam.
!!
!!         2) Calculate basic information about the linear beam ray grid.
!!
!!         3) Calculate the unit vectors associated with both the linear axis
!!            in the local target coordinate system (i.e. with origin on
!!            the center of the target line).
!!
!! ARGUMENTS
!!
!! NOTES
!!
!!***

subroutine ed_beamsInfo2DRec ()

  use Driver_interface,       ONLY : Driver_abortFlash

  use Logfile_interface,      ONLY : Logfile_stamp

  use ed_interface,           ONLY : ed_beam2DGridSetupRegular,     &
                                     ed_beam2DGridSetupStatistical, &
                                     ed_beam2DGridWeightRegular

  use EnergyDeposition_data,  ONLY : ed_beams,                 &
                                     ed_numberOfBeams,         &
                                     ed_pulseNumberOfSections, &
                                     ed_pulses

  implicit none

#include "EnergyDeposition.h"
#include "Flash.h"
#include "constants.h"

  character (len = BEAM_STRING_LENGTH) :: functionType
  character (len = BEAM_STRING_LENGTH) :: gridType

  logical :: seedIncrement, seedInitialize

  integer :: beam
  integer :: lastSection
  integer :: nRays
  integer :: nTics
  integer :: pulseNumber
  integer :: seed, seedMaximum, seedStepping

  real    :: beamLength
  real    :: beginTime, endTime
  real    :: bsizeInv
  real    :: bx, by
  real    :: delta
  real    :: firstTic
  real    :: gaussianCenter
  real    :: gaussianExponent
  real    :: gaussianRadius
  real    :: gridWeight
  real    :: lensSemiAxis
  real    :: lensX, lensY
  real    :: target2LensMagnify
  real    :: targetSemiAxis
  real    :: targetX, targetY
  real    :: ux, uy
  real    :: x,y
!
!
!     ...Loop over all beams.
!
!
  do beam = 1, ed_numberOfBeams

     functionType     = ed_beams (beam) % crossSectionFunctionType
     gaussianCenter   = ed_beams (beam) % gaussianCenterMajor
     gaussianExponent = ed_beams (beam) % gaussianExponent
     gaussianRadius   = ed_beams (beam) % gaussianRadiusMajor
     gridType         = ed_beams (beam) % gridType
     lensSemiAxis     = ed_beams (beam) % lensSemiAxisMajor
     lensX            = ed_beams (beam) % lensX
     lensY            = ed_beams (beam) % lensY
     nRays            = ed_beams (beam) % numberOfRays
     pulseNumber      = ed_beams (beam) % pulseNumber
     targetSemiAxis   = ed_beams (beam) % targetSemiAxisMajor
     targetX          = ed_beams (beam) % targetX
     targetY          = ed_beams (beam) % targetY
!
!
!     ...Set the dimensionality indicator of the beam.
!
!
     ed_beams (beam) % dimensionality = 2
!
!
!     ...Calculate distance from lens to target and pulse time bounds for each beam.
!
!
     x = targetX - lensX
     y = targetY - lensY

     beamLength = sqrt (x * x + y * y)

     lastSection = ed_pulseNumberOfSections (pulseNumber)
     beginTime   = ed_pulses (1,             pulseNumber) % pulseTime
     endTime     = ed_pulses (lastSection,   pulseNumber) % pulseTime

     ed_beams (beam) % distanceLens2Target = beamLength
     ed_beams (beam) % pulseStartingTime   = beginTime
     ed_beams (beam) % pulseEndingTime     = endTime
!
!
!     ...Calculate the magnification factor for target -> lens semiaxis. This magnification
!        factor will be needed to establish the beam grid data at the lens, since all the
!        beam grid data that has been stored is based on the target. 
!
!
     target2LensMagnify = lensSemiAxis  / targetSemiAxis

     ed_beams (beam) % target2LensMagnification = target2LensMagnify
!
!
!     ...Set up unit vectors in terms of the local target x,y-coordinate system (origin is midpoint
!        of target line).
!
!
     bx = -x         ! beam vector pointing from target to lens
     by = -y

     bSizeInv = 1.0 / beamLength

     ux = - by * bSizeInv
     uy =   bx * bSizeInv

     ed_beams (beam) % semiAxisUnitMajorX = ux
     ed_beams (beam) % semiAxisUnitMajorY = uy
!
!
!     ...Set up the beam's grid and store the info.
!
!
     if (gridType == 'regular1D') then

         call ed_beam2DGridSetupRegular  (targetsemiAxis,            &
                                          nRays,                     &
                                                            nTics,   &
                                                            delta,   &
                                                            firstTic )

         call ed_beam2DGridWeightRegular (lensX,                   &
                                          ux,                        &
                                          target2LensMagnify,        &
                                          targetSemiAxis,            &
                                          functionType,              &
                                          gaussianExponent,          &
                                          gaussianRadius,            &
                                          gaussianCenter,            &
                                          nTics,                     &
                                          delta,                     &
                                          firstTic,                  &
                                          nRays,                     &
                                                          gridWeight )
         ed_beams (beam) % gridnTics1stDim    = nTics
         ed_beams (beam) % gridDelta1stDim    = delta
         ed_beams (beam) % gridFirstTic1stDim = firstTic
         ed_beams (beam) % gridWeight         = gridWeight

     else if (gridType == 'statistical1D') then

         seedInitialize = .true.
         seedIncrement  = .false.

         call ed_beam2DGridSetupStatistical (seedInitialize,                &
                                             seedIncrement,                 &
                                                              seedMaximum,  &
                                                              seedStepping, &
                                                              seed          )

         ed_beams (beam) % gridSeed           = seed
         ed_beams (beam) % gridSeedMaximum    = seedMaximum
         ed_beams (beam) % gridSeedStepping   = seedStepping

     else
         call Driver_abortFlash ('[ed_beamsInfo2DRec] ERROR: unknown 2D beam grid type')
     end if

  end do
!
!
!     ...Ready!
!
!
  return
end subroutine ed_beamsInfo2DRec
