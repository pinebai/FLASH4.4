!!****if* source/diagnostics/ProtonImaging/localAPI/pi_xferDiskProtonsNew2OldFile
!!
!! NAME
!!
!!  pi_xferDiskProtonsNew2OldFile
!!
!! SYNOPSIS
!!
!!  call pi_xferDiskProtonsNew2OldFile (logical, intent (in) :: rewindOldFile)
!!
!! DESCRIPTION
!!
!!  When calling this routine, the complete content of the new disk proton file will
!!  be transferred to the last writing position of the old disk proton file. It can thus
!!  be viewed as appending of the new disk protons to the old disk protons. Both new and old
!!  files are assumed to exist. Only the master processor does the transfer, but the new
!!  end writing position of the old disk proton file must be known by all processors.
!!  Since this routine is called after all protons have been traced through the domain,
!!  we can use the disk proton buffer array to read in and write out the disk proton
!!  data.
!!
!! ARGUMENTS
!!
!!  rewindOldFile : If true, the old disk proton file will be written from the beginning
!!
!!***

subroutine pi_xferDiskProtonsNew2OldFile (rewindOldFile)

  implicit none

  logical, intent (in) :: rewindOldFile

  return
end subroutine pi_xferDiskProtonsNew2OldFile
