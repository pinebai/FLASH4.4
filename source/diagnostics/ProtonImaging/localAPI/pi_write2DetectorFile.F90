!!****if* source/diagnostics/ProtonImaging/localAPI/pi_write2DetectorFile
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

  implicit none

  integer, intent (in) :: detector
  integer, intent (in) :: recordCount
  integer, intent (in) :: bucketCount
  real,    intent (in) :: bucket (1:recordCount,1:bucketCount)

  return
end subroutine pi_write2DetectorFile
