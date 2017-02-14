!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/oldOneTimeStep/ProtonDetection/pi_setupDetectors
!!
!! NAME
!!
!!  pi_setupDetectors
!!
!! SYNOPSIS
!!
!!  call pi_setupDetectors ()
!!
!! DESCRIPTION
!!
!!  Sets up the detectors (location, properties) to be used for proton imaging.
!!
!! ARGUMENTS
!!
!!  none
!!
!! NOTES
!!
!!  All needed info is read in as runtime parameters.
!!
!!***

subroutine pi_setupDetectors ()

  use Driver_interface,            ONLY : Driver_abortFlash

  use pi_interface,                ONLY : pi_detectorsCheck,      &
                                          pi_detectorsInfo
  
  use ProtonImaging_data,          ONLY : pi_beamsAreSetup,      &
                                          pi_detectorFilesID,    &
                                          pi_detectors,          &
                                          pi_detectorsAreSetup,  &
                                          pi_degrees2rad,        &
                                          pi_gridGeometry,       &
                                          pi_microns2cm,         &
                                          pi_numberOfBeams,      &
                                          pi_numberOfDetectors

  use Logfile_interface,           ONLY : Logfile_stampMessage

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none

#include "ProtonImaging.h"
#include "Flash.h"
#include "constants.h"

  character (len = MAX_STRING_LENGTH) :: parameterString
  character (len = 1)                 :: sideTiltingAxis

  integer :: alignWRTbeamNr
  integer :: detector

  real    :: centerX, centerY, centerZ
  real    :: distance2BeamCapsule
  real    :: normalX, normalY, normalZ
  real    :: sideLength
  real    :: sideTiltingAngle
!
!
!     ...Allocate the detectors array.
!
!  
  allocate (pi_detectors (1:pi_numberOfDetectors))
!
!
!     ...Read the detector runtime parameters.
!
!
  do detector = 1, pi_numberOfDetectors

     write (parameterString,'(a,i0)') "pi_detectorCenterX_", detector
     call RuntimeParameters_get (parameterString, centerX)

     write (parameterString,'(a,i0)') "pi_detectorCenterY_", detector
     call RuntimeParameters_get (parameterString, centerY)

     write (parameterString,'(a,i0)') "pi_detectorCenterZ_", detector
     call RuntimeParameters_get (parameterString, centerZ)

     write (parameterString,'(a,i0)') "pi_detectorNormalX_", detector
     call RuntimeParameters_get (parameterString, normalX)

     write (parameterString,'(a,i0)') "pi_detectorNormalY_", detector
     call RuntimeParameters_get (parameterString, normalY)

     write (parameterString,'(a,i0)') "pi_detectorNormalZ_", detector
     call RuntimeParameters_get (parameterString, normalZ)

     write (parameterString,'(a,i0)') "pi_detectorSideLength_", detector
     call RuntimeParameters_get (parameterString, sideLength)

     write (parameterString,'(a,i0)') "pi_detectorSideTiltingAngle_", detector
     call RuntimeParameters_get (parameterString, sideTiltingAngle)

     write (parameterString,'(a,i0)') "pi_detectorSideTiltingAxis_", detector
     call RuntimeParameters_get (parameterString, sideTiltingAxis)

     write (parameterString,'(a,i0)') "pi_detectorAlignWRTbeamNr_", detector
     call RuntimeParameters_get (parameterString, alignWRTbeamNr)

     write (parameterString,'(a,i0)') "pi_detectorDist2BeamCapsule_", detector
     call RuntimeParameters_get (parameterString, distance2BeamCapsule)
!
!
!     ...Catch any bad data at this point.
!
!  
     if (sideLength <= 0.0) then
         write (*,*) ' Detector # , side length = ',detector, sideLength
         call Driver_abortFlash ("pi_setupDetectors: The side of the detector screen must be > 0!")
     end if

     if (sideTiltingAngle < 0.0) then
         write (*,*) ' Detector # , side tilting angle = ',detector, sideTiltingAngle
         call Driver_abortFlash ("pi_setupDetectors: The tilting angle of a side must be >= 0 degrees!")
     end if

     if (sideTiltingAngle >= 90.0) then
         write (*,*) ' Detector # , side tilting angle = ',detector, sideTiltingAngle
         call Driver_abortFlash ("pi_setupDetectors: The tilting angle of a side must be < 90 degrees!")
     end if

     if (      (sideTiltingAxis /= 'x') &
         .and. (sideTiltingAxis /= 'y') &
         .and. (sideTiltingAxis /= 'z') ) then
          write (*,*) ' Detector # , side tilting axis = ',detector, sideTiltingAxis
          call Driver_abortFlash ("ed_setupDetectors: invalid side tilting axis!")
     end if

     if (alignWRTbeamNr > pi_numberOfBeams) then
         write (*,*) ' Detector # , aligned wrt beam number = ',detector, alignWRTbeamNr
         call Driver_abortFlash ("pi_setupDetectors: Not enough beams specified!")
     end if

     if (alignWRTbeamNr > 0 .and. distance2BeamCapsule <= 0.0) then
         write (*,*) ' Detector # , distance to beam capsule center = ',detector, distance2BeamCapsule
         call Driver_abortFlash ("pi_setupDetectors: Distance to beam capsule center must be > 0!")
     end if

     if (alignWRTbeamNr > 0 .and. .not.pi_beamsAreSetup) then
         write (*,*) ' Detector # , aligned wrt beam number = ',detector, alignWRTbeamNr
         call Driver_abortFlash ("pi_setupDetectors: Beams are not set up!")
     end if
!
!
!     ...Perform some unit conversions.
!
!
     sideTiltingAngle = sideTiltingAngle * pi_degrees2rad        ! convert to radians
!
!
!     ...Store detector data obtained so far into appropriate places.
!
!  
     pi_detectors (detector) % alignWRTbeamNr            = alignWRTbeamNr
     pi_detectors (detector) % centerX                   = centerX                      ! in cm
     pi_detectors (detector) % centerY                   = centerY                      ! in cm
     pi_detectors (detector) % centerZ                   = centerZ                      ! in cm
     pi_detectors (detector) % distance2beamCapsule      = distance2beamCapsule         ! in cm
     pi_detectors (detector) % normalX                   = normalX                      ! in cm
     pi_detectors (detector) % normalY                   = normalY                      ! in cm
     pi_detectors (detector) % normalZ                   = normalZ                      ! in cm
     pi_detectors (detector) % sideLength                = sideLength                   ! in cm
     pi_detectors (detector) % sideTiltingAngle          = sideTiltingAngle             ! in rad
     pi_detectors (detector) % sideTiltingAxis           = sideTiltingAxis

  enddo
!
!
!     ...Extract more detector info (if needed) and check the gathered detectors data.
!
!  
  call pi_detectorsInfo  ()
  call pi_detectorsCheck ()
!
!
!   ...Allocate file ID array for further reference.
!
!
  allocate (pi_detectorFilesID (1:pi_numberOfDetectors))
!
!
!     ...Set detectors status indicator.
!
!
  pi_detectorsAreSetup = .true.
!
!
!     ...Ready!
!
!
  return
end subroutine pi_setupDetectors
