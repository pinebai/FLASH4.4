!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/ProtonDetection/pi_closeDetectorFiles
!!
!! NAME
!!
!!  pi_closeDetectorFiles
!!
!! SYNOPSIS
!!
!!  call pi_closeDetectorFiles ()
!!
!! DESCRIPTION
!!
!!  Closes all detector files that are currently open. Only the master processor can close
!!  detector files.
!!
!! ARGUMENTS
!!
!!***

subroutine pi_closeDetectorFiles ()

  use ProtonImaging_data,  ONLY : pi_detectorFilesID,   &
                                  pi_globalMe,          &
                                  pi_numberOfDetectors

  implicit none

#include "constants.h"
#include "ProtonImaging.h"
#include "Flash.h"

  logical :: fileOpen

  integer :: detector
  integer :: fileUnit
!
!
!   ...Only the master processor is allowed to close detector files.
!
!
  if (pi_globalMe /= MASTER_PE) then
      return
  end if
!
!
!     ...Loop over all detectors, check their file status and close if open.
!
!
  do detector = 1,pi_numberOfDetectors

     fileUnit = pi_detectorFilesID (detector)
     inquire (unit = fileUnit, opened = fileOpen)

     if (fileOpen) then
         close (fileUnit)
     end if

  end do
!
!
!    ...Ready!
!
!
  return
end subroutine pi_closeDetectorFiles
