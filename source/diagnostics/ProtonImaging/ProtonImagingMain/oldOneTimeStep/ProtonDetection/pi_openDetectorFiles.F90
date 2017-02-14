!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/oldOneTimeStep/ProtonDetection/pi_openDetectorFiles
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

  use ProtonImaging_data,  ONLY : pi_baseName,                  &
                                  pi_beams,                     &
                                  pi_detectorFileNameTimeStamp, &
                                  pi_detectorFilesID,           &
                                  pi_globalMe,                  &
                                  pi_numberOfBeams

  implicit none

#include "constants.h"
#include "ProtonImaging.h"
#include "Flash.h"

  real, intent (in) :: timeSimulation

  character (len = 2                ) :: charDetectorID
  character (len = 20               ) :: charTime
  character (len = MAX_STRING_LENGTH) :: fileName

  logical :: beamActive

  integer :: beam
  integer :: detector
  integer :: fileUnit
  integer :: nProtonsBeam
  integer :: ut_getFreeFileUnit

  real    :: time2Launch
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
!   ...Prepare the time label (just in case).
!
!
  write (charTime,'(es20.3)') timeSimulation
  charTime = adjustl (charTime)  ! remove empty characters on left
!
!
!     ...Loop over all beams, assemble the file name for each active beam and
!        open the file.
!
!
  do beam = 1, pi_numberOfBeams

     time2Launch  = pi_beams (beam) % time2Launch
     nProtonsBeam = pi_beams (beam) % numberOfProtons

     beamActive = (time2Launch <= timeSimulation) .and. (nProtonsBeam > 0)

     if (beamActive) then

         detector = pi_beams (beam) % detector
         write (charDetectorID,'(I2.2)') detector

         fileUnit = ut_getFreeFileUnit ()

         fileName = trim (pi_baseName)   // &
                    "ProtonDetectorFile" // &
                    charDetectorID

         if (pi_detectorFileNameTimeStamp) then
             fileName = trim (fileName) // "_" // trim (charTime)
         end if

         open (fileUnit, file = fileName, status = 'REPLACE')

         pi_detectorFilesID (detector) = fileUnit

     end if
  end do
!
!
!    ...Ready!
!
!
  return
end subroutine pi_openDetectorFiles
