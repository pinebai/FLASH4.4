!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/ProtonDetection/pi_setupDetectorFileNames
!!
!! NAME
!!
!!  pi_setupDetectorFileNames
!!
!! SYNOPSIS
!!
!!  call pi_setupDetectorFileNames (real, intent (in) :: timeSimulation)
!!
!! DESCRIPTION
!!
!!  Sets the detector file names for further use during a proton imaging step. The
!!  names of the detector files are as follows:
!!
!!              <basenm> ProtonDetectorFile <detectorID> _ <timeSimulation>
!!
!!  where <basenm> is the simulation base name, <detectorID> contains the detector number
!!  and (optionally) <timeSimulation> is the current simulation time.
!!
!! ARGUMENTS
!!
!!  timeSimulation  : current simulation time
!!
!! NOTES
!!          
!!  1) As many detector file names are established as there were number of detectors
!!     specified. The detector number range is currently limited to 01-99. If more
!!     detectors were specified, the routine aborts.
!!
!!  2) A runtime parameter exists which enables suppression of the _<timeSimulation> part
!!     of the detector file names.
!!
!!  3) Although only the master processor will ever need these names, all processors
!!     get a copy of the detector names.
!!
!!***

subroutine pi_setupDetectorFileNames (timeSimulation)

  use Driver_interface,    ONLY : Driver_abortFlash

  use ProtonImaging_data,  ONLY : pi_baseName,                  &
                                  pi_detectorFileNameTimeStamp, &
                                  pi_detectorFilesName,         &
                                  pi_globalMe,                  &
                                  pi_numberOfDetectors

  implicit none

#include "constants.h"
#include "ProtonImaging.h"
#include "Flash.h"

  real, intent (in) :: timeSimulation

  character (len = 2                ) :: charDetectorNumber
  character (len = 20               ) :: charTime
  character (len = MAX_STRING_LENGTH) :: fileName

  integer :: detector
!
!
!   ...Abort, if number of detectors > 99.
!
!
  if (pi_numberOfDetectors > 99) then
      call Driver_abortFlash ("pi_setupDetectorFileNames: number of detectors > 99!")
  end if
!
!
!   ...Prepare the time label (just in case).
!
!
  write (charTime,'(es20.3)') timeSimulation
  charTime = adjustl (charTime)  ! remove empty characters on left
!
!
!     ...Loop over all detectors and assemble the file name for each detector.
!
!
  do detector = 1, pi_numberOfDetectors

     write (charDetectorNumber,'(I2.2)') detector

     fileName = trim (pi_baseName)   // "ProtonDetectorFile" // charDetectorNumber

     if (pi_detectorFileNameTimeStamp) then
         fileName = trim (fileName) // "_" // trim (charTime)
     end if

     pi_detectorFilesName (detector) = fileName

  end do
!
!
!    ...Ready!
!
!
  return
end subroutine pi_setupDetectorFileNames
