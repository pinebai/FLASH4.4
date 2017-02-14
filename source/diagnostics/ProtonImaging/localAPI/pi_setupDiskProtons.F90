!!****if* source/diagnostics/ProtonImaging/localAPI/pi_setupDiskProtons
!!
!! NAME
!!
!!  pi_setupDiskProtons
!!
!! SYNOPSIS
!!
!!  call pi_setupDiskProtons ()
!!
!! DESCRIPTION
!!
!!  Sets up the disk protons, which means to allocate the needed disk proton array.
!!  The disk proton array will be used to store protons that remain in the domain
!!  during a time step and will be written to disk for consideration during the
!!  following time step. The old and new disk proton file names are set here. The
!!  disk proton file names are currently as follows:
!!
!!                          <basenm> ProtonImagingDiskProtonsNew
!!                          <basenm> ProtonImagingDiskProtonsOld
!!
!!  where <basenm> is the simulation base name. The routine also sets up and commits
!!  a mpi type structure corresponding to the disk protons, which will be used when
!!  sending disk protons between processors.
!!
!! ARGUMENTS
!!
!! NOTES
!!
!!***

subroutine pi_setupDiskProtons ()

  implicit none

  return
end subroutine pi_setupDiskProtons
