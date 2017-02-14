!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserBeams/ed_initializeBeams
!!
!! NAME
!!
!!  ed_initializeBeams
!!
!! SYNOPSIS
!!
!!  call ed_initializeBeams ()
!!
!! DESCRIPTION
!!
!!  Initializes the beams.
!!
!! ARGUMENTS
!!
!! NOTES
!!
!!***

subroutine ed_initializeBeams ()

  use EnergyDeposition_data,   ONLY : ed_beams,         &
                                      ed_notSetInteger, &
                                      ed_notSetReal,    &
                                      ed_numberOfBeams
  implicit none

  integer :: beam
!
!
!     ...Do the initialization.
!
!
  do beam = 1, ed_numberOfBeams

     ed_beams (beam) % crossSectionFunctionType  = ' '
     ed_beams (beam) % dimensionality            = ed_notSetInteger
     ed_beams (beam) % distanceLens2Target       = ed_notSetReal
     ed_beams (beam) % frequency                 = ed_notSetReal
     ed_beams (beam) % gaussianCenterMajor       = ed_notSetReal
     ed_beams (beam) % gaussianCenterMinor       = ed_notSetReal
     ed_beams (beam) % gaussianExponent          = ed_notSetReal
     ed_beams (beam) % gaussianRadiusMajor       = ed_notSetReal
     ed_beams (beam) % gaussianRadiusMinor       = ed_notSetReal
     ed_beams (beam) % gridDelta1stDim           = ed_notSetReal
     ed_beams (beam) % gridDelta2ndDim           = ed_notSetReal
     ed_beams (beam) % gridFirstTic1stDim        = ed_notSetReal
     ed_beams (beam) % gridFirstTic2ndDim        = ed_notSetReal
     ed_beams (beam) % gridnTics1stDim           = ed_notSetInteger
     ed_beams (beam) % gridnTics2ndDim           = ed_notSetInteger
     ed_beams (beam) % gridSeed                  = ed_notSetInteger
     ed_beams (beam) % gridSeedMaximum           = ed_notSetInteger
     ed_beams (beam) % gridSeedStepping          = ed_notSetInteger
     ed_beams (beam) % gridType                  = ' '
     ed_beams (beam) % gridWeight                = ed_notSetReal
     ed_beams (beam) % ignoreBoundaryCondition   = .false.
     ed_beams (beam) % initialRaySpeed           = ed_notSetReal
     ed_beams (beam) % lensSemiAxisMajor         = ed_notSetReal
     ed_beams (beam) % lensSemiAxisMinor         = ed_notSetReal
     ed_beams (beam) % lensX                     = ed_notSetReal
     ed_beams (beam) % lensY                     = ed_notSetReal
     ed_beams (beam) % lensZ                     = ed_notSetReal
     ed_beams (beam) % numberOfRays              = ed_notSetInteger
     ed_beams (beam) % pulseNumber               = ed_notSetInteger
     ed_beams (beam) % pulseStartingTime         = ed_notSetReal
     ed_beams (beam) % pulseEndingTime           = ed_notSetReal
     ed_beams (beam) % semiAxisMajorTorsionAngle = ed_notSetReal
     ed_beams (beam) % semiAxisMajorTorsionAxis  = ' '
     ed_beams (beam) % semiAxisUnitMajorX        = ed_notSetReal
     ed_beams (beam) % semiAxisUnitMajorY        = ed_notSetReal
     ed_beams (beam) % semiAxisUnitMajorZ        = ed_notSetReal
     ed_beams (beam) % semiAxisUnitMinorX        = ed_notSetReal
     ed_beams (beam) % semiAxisUnitMinorY        = ed_notSetReal
     ed_beams (beam) % semiAxisUnitMinorZ        = ed_notSetReal
     ed_beams (beam) % target2LensMagnification  = ed_notSetReal
     ed_beams (beam) % targetSemiAxisMajor       = ed_notSetReal
     ed_beams (beam) % targetSemiAxisMinor       = ed_notSetReal
     ed_beams (beam) % targetX                   = ed_notSetReal
     ed_beams (beam) % targetY                   = ed_notSetReal
     ed_beams (beam) % targetZ                   = ed_notSetReal
     ed_beams (beam) % wavelength                = ed_notSetReal

  end do
!
!
!     ...Ready!
!
!
  return
end subroutine ed_initializeBeams
