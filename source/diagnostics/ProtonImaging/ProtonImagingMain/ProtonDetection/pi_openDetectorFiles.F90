!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/ProtonDetection/pi_openDetectorFiles
!!
!! NAME
!!
!!  pi_openDetectorFiles
!!
!! SYNOPSIS
!!
!!  call pi_openDetectorFiles (real, intent (in) :: timeSimulation)
!!
!! DESCRIPTION
!!
!!  Opens the detector files for recording screen protons for active proton beams only.
!!  To each detector there corresponds one file, characterized by the detector number and
!!  the current simulation time. Only the master processor opens the detector files. The
!!  names of the detector files is as follows:
!!
!!              <basenm> ProtonDetectorFile <detectorID> _ <timeSimulation>
!!
!!  where <basenm> is the simulation base name, <detectorID> contains the detector number
!!  and <timeSimulation> is the current simulation time.
!!
!! ARGUMENTS
!!
!!  timeSimulation  : current simulation time
!!
!! NOTES
!!          
!!  1) Only detector files are opened for active beams. The detector number range is
!!     currently limited to 01-99. If this routine is called, then there is at least one
!!     active proton beam.
!!
!!  2) A runtime parameter exists which enables suppression of the _<timeSimulation> part
!!     of the detector file names.
!!
!!***

subroutine pi_openDetectorFiles (timeSimulation)

  use Driver_interface,    ONLY : Driver_abortFlash

  use ProtonImaging_data,  ONLY : pi_detectorFilesID,           &
                                  pi_detectorFilesName,         &
                                  pi_globalMe,                  &
                                  pi_numberOfDetectors

  implicit none

#include "constants.h"
#include "ProtonImaging.h"
#include "Flash.h"

  real, intent (in) :: timeSimulation

  character (len = MAX_STRING_LENGTH) :: fileName

  logical :: fileExists
  logical :: fileOpen

  integer :: detector
  integer :: fileUnit
  integer :: ut_getFreeFileUnit
!
!
!   ...Open the files only on the master processor.
!
!
  if (pi_globalMe /= MASTER_PE) then
      return
  end if
!
!
!     ...Looop over all detectors and open the files.
!
!
  do detector = 1, pi_numberOfDetectors

      fileName = pi_detectorFilesName (detector)

      inquire (file = fileName, exist = fileExists, opened = fileOpen)

      if (fileOpen) then
          call Driver_abortFlash ("pi_openDetectorFile: detector file already open!")
      end if

      fileUnit = ut_getFreeFileUnit ()

      if (fileExists) then
          open (fileUnit, file = fileName, status = 'OLD', position = 'APPEND')
      else
          open (fileUnit, file = fileName, status = 'NEW')
      end if

      pi_detectorFilesID (detector) = fileUnit

  end do
!
!
!    ...Ready!
!
!
  return
end subroutine pi_openDetectorFiles
