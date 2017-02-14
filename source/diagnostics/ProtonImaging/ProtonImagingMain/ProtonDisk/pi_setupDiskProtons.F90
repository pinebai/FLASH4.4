!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/ProtonDisk/pi_setupDiskProtons
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

  use ProtonImaging_data,  ONLY : pi_baseName,                  &
                                  pi_diskProtonNewFileName,     &
                                  pi_diskProtonNewFileNbatches, &
                                  pi_diskProtonOldFileName,     &
                                  pi_diskProtonOldFileNbatches, &
                                  pi_diskProtons,               &
                                  pi_maxProtonCount,            &
                                  pi_mpiDiskProtonType

  implicit none

#include "ProtonImaging.h"
 include "Flash_mpi.h"

  integer :: error
!
!
!     ...Allocate the disk proton array and set the disk proton file names
!        and file IDs. Initialize the disk proton batch count on both files.
!
!
  allocate (pi_diskProtons (1:PROTON_ATTRCOUNT,1:pi_maxProtonCount))

  pi_diskProtonNewFileName = trim (pi_baseName) // "ProtonImagingDiskProtonsNew"
  pi_diskProtonOldFileName = trim (pi_baseName) // "ProtonImagingDiskProtonsOld"

  pi_diskProtonNewFileNbatches = 0
  pi_diskProtonOldFileNbatches = 0
!
!
!     ...Create and commit the disk proton mpi type structure.
!
!
  call MPI_Type_Contiguous (PROTON_ATTRCOUNT,                 &
                            FLASH_REAL,                       &
                                        pi_mpiDiskProtonType, &
                                        error                 )

  call MPI_Type_Commit (pi_mpiDiskProtonType, error)
!
!
!     ...Ready!
!
!
  return
end subroutine pi_setupDiskProtons
