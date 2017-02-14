!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/ProtonDisk/pi_xferDiskProtonsNew2OldFile
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
!!  number of disk proton batches of the old disk proton file must be known by all processors.
!!  Since this routine is called after all protons have been traced through the domain,
!!  we can use the disk proton buffer array to read in and write out the disk proton
!!  data.
!!
!! ARGUMENTS
!!
!!  rewindOldFile : If true, the old disk proton file is written from the beginning
!!
!!***

subroutine pi_xferDiskProtonsNew2OldFile (rewindOldFile)

  use Driver_interface,    ONLY : Driver_abortFlash

  use ProtonImaging_data,  ONLY : pi_diskProtonNewFileID,       &
                                  pi_diskProtonNewFileNbatches, &
                                  pi_diskProtonOldFileID,       &
                                  pi_diskProtonOldFileNbatches, &
                                  pi_diskProtons,               &
                                  pi_globalComm,                &
                                  pi_globalMe,                  &
                                  pi_maxProtonCount,            &
                                  pi_monitorFileUnit
                                  

  implicit none

#include "constants.h"
#include "ProtonImaging.h"
#include "Flash.h"

  include "Flash_mpi.h"

  logical, intent (in) :: rewindOldFile

  logical :: fileExists

  integer :: batch
  integer :: error
  integer :: numberOfProtons
  integer :: xferNumberOfProtons
!
!
!   ...Check first for existence of both files (master only).
!
!
  if (pi_globalMe == MASTER_PE) then

      inquire (unit = pi_diskProtonOldFileID, exist = fileExists)

      if (.not.fileExists) then
           call Driver_abortFlash ("pi_xferDiskProtonsNew2OldFile: Old disk proton file does not exist!")
      end if

      inquire (unit = pi_diskProtonNewFileID, exist = fileExists)

      if (.not.fileExists) then
           call Driver_abortFlash ("pi_xferDiskProtonsNew2OldFile: New disk proton file does not exist!")
      end if
!
!
!   ...Make the transfer by reading/writing batches of disk protons (master only).
!      Update the old disk proton file counter locally on the master.
!
!
      if (rewindOldFile) then
          rewind (pi_diskProtonOldFileID)
          pi_diskProtonOldFileNbatches = 0
      end if

      rewind (pi_diskProtonNewFileID)

      xferNumberOfProtons = 0

      if (pi_diskProtonNewFileNbatches > 0) then
          do batch = 1,pi_diskProtonNewFileNbatches
             read  (pi_diskProtonNewFileID) numberOfProtons
             if (numberOfProtons > 0) then
                 read  (pi_diskProtonNewFileID) pi_diskProtons (1:PROTON_ATTRCOUNT,1:numberOfProtons)
                 write (pi_diskProtonOldFileID)  numberOfProtons
                 write (pi_diskProtonOldFileID) pi_diskProtons (1:PROTON_ATTRCOUNT,1:numberOfProtons)
                 xferNumberOfProtons = xferNumberOfProtons + numberOfProtons
                 pi_diskProtonOldFileNbatches = pi_diskProtonOldFileNbatches + 1
             end if
          end do
      end if

      if (xferNumberOfProtons > 0) then
          write (pi_monitorFileUnit,'(a,i8,a,i8,a)') ' xferred ', xferNumberOfProtons, &
                 ' New    Disk   Protons -> Old Disk (total # proton batches = ',pi_diskProtonOldFileNbatches,')'
      end if

  end if
!
!
!   ...Share the updated batch counter of the old disk proton file with all processors.
!
!
  call MPI_Bcast (pi_diskProtonOldFileNbatches, &
                  1,                            &
                  MPI_INTEGER,                  &
                  MASTER_PE,                    &
                  pi_globalComm,                &
                  error                         )
!
!
!    ...Ready!
!
!
  return
end subroutine pi_xferDiskProtonsNew2OldFile
