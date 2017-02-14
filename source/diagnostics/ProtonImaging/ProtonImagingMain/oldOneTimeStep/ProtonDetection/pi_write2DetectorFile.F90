!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/oldOneTimeStep/ProtonDetection/pi_write2DetectorFile
!!
!! NAME
!!
!!  pi_write2DetectorFile
!!
!! SYNOPSIS
!!
!!  call pi_write2DetectorFile (integer, intent (in) :: detector,
!!                              integer, intent (in) :: recordCount,
!!                              integer, intent (in) :: bucketCount,
!!                              real,    intent (in) :: bucket (:,:))
!!
!! DESCRIPTION
!!
!!  Writes a bucket of screen protons to the appropriate detector file on disk. The file is
!!  already open and the screen protons data is appended. Only the master processor writes
!!  out the bucket.
!!
!! ARGUMENTS
!!
!!  detector      : The detector number
!!  recordCount   : Number of records in bucket (up to 6: screen [x,y] pairs + 4 diagnostics)
!!  bucketCount   : Number of elements in bucket
!!  bucket        : The bucket containing all the elements.
!!
!!***

subroutine pi_write2DetectorFile (detector, recordCount, bucketCount, bucket)

  use Driver_interface,    ONLY : Driver_abortFlash

  use ProtonImaging_data,  ONLY : pi_detectorFilesID,       &
                                  pi_detectorLNwriteFormat, &
                                  pi_globalMe

  implicit none

#include "constants.h"
#include "ProtonImaging.h"
#include "Flash.h"

  integer, intent (in) :: detector
  integer, intent (in) :: recordCount
  integer, intent (in) :: bucketCount
  real,    intent (in) :: bucket (1:recordCount,1:bucketCount)

  integer :: fileUnit
!
!
!   ...Write out the buckets only on the master processor.
!
!
  if (pi_globalMe == MASTER_PE) then
      fileUnit = pi_detectorFilesID (detector)
      write (fileUnit, pi_detectorLNwriteFormat) (bucket (1:recordCount,1:bucketCount))
  end if
!
!
!    ...Ready!
!
!
  return
end subroutine pi_write2DetectorFile
