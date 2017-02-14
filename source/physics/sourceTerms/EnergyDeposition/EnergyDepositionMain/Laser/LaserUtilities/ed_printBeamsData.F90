!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserUtilities/ed_printBeamsData
!!
!! NAME
!!
!!  ed_printBeamsData
!!
!! SYNOPSIS
!!
!!  call ed_printBeamsData ()
!!
!! DESCRIPTION
!!
!!  Utility routine, which prints detailed info regarding the generated data
!!  for all beams to a text file. The information is written out to a file named
!!  <basenm>LaserBeamsPrint.txt, where <basenm> is the runtime parameter for
!!  output file names.
!!
!!***

subroutine ed_printBeamsData ()

  use Driver_interface,       ONLY : Driver_abortFlash

  use EnergyDeposition_data,  ONLY : ed_baseName, &
                                     ed_beams,    &
                                     ed_globalMe, &
                                     ed_numberOfBeams

  implicit none
   
#include "constants.h"
#include "EnergyDeposition.h"
#include "Flash.h"

  character (len =  MAX_STRING_LENGTH) :: fileName
  character (len = BEAM_STRING_LENGTH) :: functionType
  character (len = BEAM_STRING_LENGTH) :: gridType
  character (len = 1)                  :: torsionAxis

  logical :: ignoreBnd

  integer :: beam
  integer :: dimensions
  integer :: fileUnit
  integer :: nRays
  integer :: nTics1D, nTics2D
  integer :: pulseNumber
  integer :: seed, seedMax, seedStep
  integer :: ut_getFreeFileUnit

  real    :: axisMajL, axisMinL
  real    :: axisMajT, axisMinT
  real    :: gaussCenMaj, gaussCenMin
  real    :: gaussExp
  real    :: gaussRadMaj, gaussRadMin
  real    :: delta1D, delta2D
  real    :: distL2T
  real    :: firstTic1D, firstTic2D
  real    :: frequency
  real    :: gridWeight
  real    :: lensX, lensY, lensZ
  real    :: magnifyT2L
  real    :: pulseStart, pulseEnd
  real    :: raySpeed
  real    :: targetX, targetY, targetZ
  real    :: torsionAngle
  real    :: unit1X, unit1Y, unit1Z
  real    :: unit2X, unit2Y, unit2Z
  real    :: wavelength
!
!
!   ...Do the printout only on the master processor.
!
!
  if (ed_globalMe /= MASTER_PE) then
      return
  end if
!
!
!   ...Open the printout file.
!
!
  fileUnit = ut_getFreeFileUnit ()
  fileName = trim (ed_baseName) // "LaserBeamsPrint.txt"

  open (fileUnit, file = fileName)
!
!
!     ...Loop over all beams.
!
!
  do beam = 1, ed_numberOfBeams

     functionType = ed_beams (beam) % crossSectionFunctionType
     dimensions   = ed_beams (beam) % dimensionality
     distL2T      = ed_beams (beam) % distanceLens2Target
     frequency    = ed_beams (beam) % frequency
     gaussCenMaj  = ed_beams (beam) % gaussianCenterMajor
     gaussCenMin  = ed_beams (beam) % gaussianCenterMinor
     gaussExp     = ed_beams (beam) % gaussianExponent
     gaussRadMaj  = ed_beams (beam) % gaussianRadiusMajor
     gaussRadMin  = ed_beams (beam) % gaussianRadiusMinor
     delta1D      = ed_beams (beam) % gridDelta1stDim
     delta2D      = ed_beams (beam) % gridDelta2ndDim
     firstTic1D   = ed_beams (beam) % gridFirstTic1stDim
     firstTic2D   = ed_beams (beam) % gridFirstTic2ndDim
     nTics1D      = ed_beams (beam) % gridnTics1stDim
     nTics2D      = ed_beams (beam) % gridnTics2ndDim
     seed         = ed_beams (beam) % gridSeed
     seedMax      = ed_beams (beam) % gridSeedMaximum
     seedStep     = ed_beams (beam) % gridSeedStepping
     gridType     = ed_beams (beam) % gridType
     gridWeight   = ed_beams (beam) % gridWeight
     ignoreBnd    = ed_beams (beam) % ignoreBoundaryCondition
     raySpeed     = ed_beams (beam) % initialRaySpeed
     axisMajL     = ed_beams (beam) % lensSemiAxisMajor
     axisMinL     = ed_beams (beam) % lensSemiAxisMinor
     lensX        = ed_beams (beam) % lensX
     lensY        = ed_beams (beam) % lensY
     lensZ        = ed_beams (beam) % lensZ
     nRays        = ed_beams (beam) % numberOfRays
     pulseNumber  = ed_beams (beam) % pulseNumber
     pulseStart   = ed_beams (beam) % pulseStartingTime
     pulseEnd     = ed_beams (beam) % pulseEndingTime
     torsionAngle = ed_beams (beam) % semiAxisMajorTorsionAngle
     torsionAxis  = ed_beams (beam) % semiAxisMajorTorsionAxis
     unit1X       = ed_beams (beam) % semiAxisUnitMajorX
     unit1Y       = ed_beams (beam) % semiAxisUnitMajorY
     unit1Z       = ed_beams (beam) % semiAxisUnitMajorZ
     unit2X       = ed_beams (beam) % semiAxisUnitMinorX
     unit2Y       = ed_beams (beam) % semiAxisUnitMinorY
     unit2Z       = ed_beams (beam) % semiAxisUnitMinorZ
     magnifyT2L   = ed_beams (beam) % target2LensMagnification
     axisMajT     = ed_beams (beam) % targetSemiAxisMajor
     axisMinT     = ed_beams (beam) % targetSemiAxisMinor
     targetX      = ed_beams (beam) % targetX
     targetY      = ed_beams (beam) % targetY
     targetZ      = ed_beams (beam) % targetZ
     wavelength   = ed_beams (beam) % wavelength

     write (fileUnit,'(/)'        )
     write (fileUnit,'(a,i2)'     ) "                LASER BEAM NR ",beam
     write (fileUnit,'(/)'        )
     write (fileUnit,'(a,es20.12)') " Distance lens center --> target center         = ", distL2T
     write (fileUnit,'(a,es20.12)') " Frequency of laser beam                        = ", frequency
     write (fileUnit,'(a,es20.12)') " Wavelength of laser beam                       = ", wavelength
     write (fileUnit,'(a,es20.12)') " Initial ray speed (in units of light speed)    = ", raySpeed
     write (fileUnit,'(a,i6)'     ) " Pulse identification number                    = ", pulseNumber
     write (fileUnit,'(a,es20.12)') " Pulse overall starting time                    = ", pulseStart
     write (fileUnit,'(a,es20.12)') " Pulse overall ending time                      = ", pulseEnd
     write (fileUnit,'(a,i20)'    ) " Number of rays in beam                         = ", nRays
     write (fileUnit,'(a,L1)'     ) " The domain boundary conditions are ignored     = ", ignoreBnd

     select case (dimensions)

     case (1)

     write (fileUnit,'(a,es20.12)') " Lens global x-coordinate                       = ", lensX
     write (fileUnit,'(a,es20.12)') " Target global x-coordinate                     = ", targetX
     write (fileUnit,'(a)'        ) " ------------------------------------------------ "

     case (2)

     write (fileUnit,'(a,a)'      ) " Cross section weight function type             = ", functionType
     write (fileUnit,'(a,a)'      ) " Grid type                                      = ", gridType
     write (fileUnit,'(a,es20.12)') " Grid weight                                    = ", gridWeight
     write (fileUnit,'(a,es20.12)') " Gaussian exponent                              = ", gaussExp
     write (fileUnit,'(a,es20.12)') " Gaussian center location along semiaxis        = ", gaussCenMaj
     write (fileUnit,'(a,es20.12)') " Gaussian radius along semiaxis                 = ", gaussRadMaj
     write (fileUnit,'(a,i20)'    ) " Number of tics on linear grid                  = ", nTics1D
     write (fileUnit,'(a,i20)'    ) " Current seed value for random grid             = ", seed
     write (fileUnit,'(a,i20)'    ) " Maximum seed value for random grid             = ", seedMax
     write (fileUnit,'(a,i20)'    ) " Stepping seed value for random grid            = ", seedStep
     write (fileUnit,'(a,es20.12)') " Tic separation on linear grid                  = ", delta1D
     write (fileUnit,'(a,es20.12)') " First tic position on linear grid              = ", firstTic1D
     write (fileUnit,'(a,es20.12)') " Lens linear semiaxis length                    = ", axisMajL
     write (fileUnit,'(a,es20.12)') " Lens global x-coordinate                       = ", lensX
     write (fileUnit,'(a,es20.12)') " Lens global y-coordinate                       = ", lensY
     write (fileUnit,'(a,es20.12)') " Linear semiaxis unit vector x-coordinate       = ", unit1X
     write (fileUnit,'(a,es20.12)') " Linear semiaxis unit vector y-coordinate       = ", unit1Y
     write (fileUnit,'(a,es20.12)') " Target linear semiaxis length                  = ", axisMajT
     write (fileUnit,'(a,es20.12)') " Target global x-coordinate                     = ", targetX
     write (fileUnit,'(a,es20.12)') " Target global y-coordinate                     = ", targetY
     write (fileUnit,'(a,es20.12)') " Target -> lens magnification factor            = ", magnifyT2L
     write (fileUnit,'(a)'        ) " ------------------------------------------------ "

     case (3)

     write (fileUnit,'(a,a)'      ) " Cross section weight function type             = ", functionType
     write (fileUnit,'(a,a)'      ) " Grid type                                      = ", gridType
     write (fileUnit,'(a,es20.12)') " Grid weight                                    = ", gridWeight
     write (fileUnit,'(a,es20.12)') " Gaussian exponent                              = ", gaussExp
     write (fileUnit,'(a,es20.12)') " Gaussian center location along major semiaxis  = ", gaussCenMaj
     write (fileUnit,'(a,es20.12)') " Gaussian center location along minor semiaxis  = ", gaussCenMin
     write (fileUnit,'(a,es20.12)') " Gaussian radius along major semiaxis           = ", gaussRadMaj
     write (fileUnit,'(a,es20.12)') " Gaussian radius along minor semiaxis           = ", gaussRadMin
     write (fileUnit,'(a,i20)'    ) " Number of tics on grids 1st dimension          = ", nTics1D
     write (fileUnit,'(a,i20)'    ) " Number of tics on grids 2nd dimension          = ", nTics2D
     write (fileUnit,'(a,i20)'    ) " Current seed value for random grid             = ", seed
     write (fileUnit,'(a,i20)'    ) " Maximum seed value for random grid             = ", seedMax
     write (fileUnit,'(a,i20)'    ) " Stepping seed value for random grid            = ", seedStep
     write (fileUnit,'(a,es20.12)') " Tic separation on grid 1st dimension           = ", delta1D
     write (fileUnit,'(a,es20.12)') " Tic separation on grid 2nd dimension           = ", delta2D
     write (fileUnit,'(a,es20.12)') " First tic position on grid 1st dimension       = ", firstTic1D
     write (fileUnit,'(a,es20.12)') " First tic position on grid 2nd dimension       = ", firstTic2D
     write (fileUnit,'(a,es20.12)') " Lens major semiaxis length                     = ", axisMajL
     write (fileUnit,'(a,es20.12)') " Lens minor semiaxis length                     = ", axisMinL
     write (fileUnit,'(a,es20.12)') " Lens global x-coordinate                       = ", lensX
     write (fileUnit,'(a,es20.12)') " Lens global y-coordinate                       = ", lensY
     write (fileUnit,'(a,es20.12)') " Lens global z-coordinate                       = ", lensZ
     write (fileUnit,'(a,es20.12)') " Major semiaxis torsion angle (rad)             = ", torsionAngle
     write (fileUnit,'(a,a1)'     ) " Major semiaxis torsion axis                    = ", torsionAxis
     write (fileUnit,'(a,es20.12)') " Major semiaxis unit vector x-coordinate        = ", unit1X
     write (fileUnit,'(a,es20.12)') " Major semiaxis unit vector y-coordinate        = ", unit1Y
     write (fileUnit,'(a,es20.12)') " Major semiaxis unit vector z-coordinate        = ", unit1Z
     write (fileUnit,'(a,es20.12)') " Minor semiaxis unit vector x-coordinate        = ", unit2X
     write (fileUnit,'(a,es20.12)') " Minor semiaxis unit vector y-coordinate        = ", unit2Y
     write (fileUnit,'(a,es20.12)') " Minor semiaxis unit vector z-coordinate        = ", unit2Z
     write (fileUnit,'(a,es20.12)') " Target major semiaxis length                   = ", axisMajT
     write (fileUnit,'(a,es20.12)') " Target minor semiaxis length                   = ", axisMinT
     write (fileUnit,'(a,es20.12)') " Target global x-coordinate                     = ", targetX
     write (fileUnit,'(a,es20.12)') " Target global y-coordinate                     = ", targetY
     write (fileUnit,'(a,es20.12)') " Target global z-coordinate                     = ", targetZ
     write (fileUnit,'(a,es20.12)') " Target -> lens magnification factor            = ", magnifyT2L
     write (fileUnit,'(a)'        ) " ------------------------------------------------ "

     case default

     call Driver_abortFlash ("ed_printBeamsData: Bad dimensionality of beam (programmers fault)!")

     end select

  end do
!
!
!   ...Close the printout file.
!
!
  close (fileUnit)
!
!
!    ...Ready!
!
!  
  return
end subroutine ed_printBeamsData
