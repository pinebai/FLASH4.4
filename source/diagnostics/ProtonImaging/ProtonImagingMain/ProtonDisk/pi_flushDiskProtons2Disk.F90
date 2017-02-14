!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/ProtonDisk/pi_flushDiskProtons2Disk
!!
!! NAME
!!
!!  pi_flushDiskProtons2Disk
!!
!! SYNOPSIS
!!
!!  call pi_flushDiskProtons2Disk ()
!!
!! DESCRIPTION
!!
!!  When calling this routine, the disk protons accumulated by each processor are all
!!  send to the master processor and processed for writing out to disk. The following
!!  steps are performed:
!!
!!     1) Gather at the master processor the info of how many disk protons each
!!        processor currently has. Calculate the offsets to prepare for storage
!!        of the disk protons on the master processor.
!!
!!     2) Gather all the disk protons on the master processor. The disk protons
!!        on the master remain in position, the others are appended according to the
!!        offsets calculated previously.
!!
!!     3) On the master processor, write out the complete set of disk protons to disk
!!        either by appending them to an already existing disk proton file or by
!!        creating a new disk proton file.
!!
!!     4) Once all disk protons have been written to file, reset the disk proton counter
!!        to 0 on all processors.
!!
!! ARGUMENTS
!!
!! NOTES
!!
!!***

subroutine pi_flushDiskProtons2Disk ()

  use ProtonImaging_data,  ONLY : pi_diskProtonCount,           &
                                  pi_diskProtonCountOffsets,    &
                                  pi_diskProtonCountProcs,      &
                                  pi_diskProtonNewFileID,       &
                                  pi_diskProtonNewFileNbatches, &
                                  pi_diskProtons,               &
                                  pi_globalComm,                &
                                  pi_globalDiskProtonCount,     &
                                  pi_globalMe,                  &
                                  pi_globalNumProcs,            &
                                  pi_maxProtonCount,            &
                                  pi_monitorFileUnit,           &
                                  pi_mpiDiskProtonType
  
  use Driver_interface,    ONLY : Driver_abortFlash

  implicit none

#include "Flash.h"
#include "constants.h"
#include "ProtonImaging.h"

  include "Flash_mpi.h"

  logical :: fileExists

  integer :: error
  integer :: offset
  integer :: rank
!
!
!     ...Gather the number of disk protons from all processors at the master processor.
!
!
  call MPI_Gather (pi_diskProtonCount,         &       ! this is being sent from rank i process
                   1,                          &       ! how many elements ?
                   MPI_INTEGER,                &       ! type of elements sent
                   pi_diskProtonCountProcs,    &       ! stored into i-th array element on master only
                   1,                          &       ! number of elements stored (master)
                   MPI_INTEGER,                &       ! type of elements stored (master)
                   MASTER_PE,                  &       ! rank of master
                   pi_globalComm,              &       ! communicator handle
                   error                       )       ! error handle
!
!
!     ...Calculate the offsets on the master only. The offset for the master is set
!        as zero, as we do not want to move disk protons on the master. Check, if
!        the total number of expected disk protons does not exceed the maximum
!        declared storage.
!
!        Next gather all disk protons from all processors on the master processor.
!        The master at this point has already the offset info for each processor.
!        Do not move disk protons on the master.
!
!
  if (pi_globalMe == MASTER_PE) then

      pi_diskProtonCountOffsets (MASTER_PE) = 0

      offset = pi_diskProtonCountProcs (MASTER_PE)

      do rank = 0,MASTER_PE-1
         pi_diskProtonCountOffsets (rank) = offset
         offset = offset + pi_diskProtonCountProcs (rank)
      end do

      do rank = MASTER_PE+1,pi_globalNumProcs-1
         pi_diskProtonCountOffsets (rank) = offset
         offset = offset + pi_diskProtonCountProcs (rank)
      end do

      pi_globalDiskProtonCount = offset   ! this is the new final (total) disk proton count on the master

      if (pi_globalDiskProtonCount > pi_maxProtonCount) then
          call Driver_abortFlash ('[pi_flushDiskProtons2Disk] ERROR: global # of disk protons > storage size')
      end if

      call MPI_Gatherv (MPI_IN_PLACE,               &       ! do not move disk protons on master
                        pi_diskProtonCount,         &       ! irrelevant
                        pi_mpiDiskProtonType,       &       ! irrelevant
                        pi_diskProtons,             &       ! stored into array on master only
                        pi_diskProtonCountProcs,    &       ! number of elements stored (master)
                        pi_diskProtonCountOffsets,  &       ! where elements are stored (master)
                        pi_mpiDiskProtonType,       &       ! type of elements stored (master)
                        MASTER_PE,                  &       ! rank of master
                        pi_globalComm,              &       ! communicator handle
                        error                       )       ! error handle

      if (pi_globalDiskProtonCount > 0) then
          write (pi_monitorFileUnit,'(a,i8,a)') ' stored  ',pi_globalDiskProtonCount,' Domain Disk   Protons'
      end if

  else

      call MPI_Gatherv (pi_diskProtons,             &       ! this is being sent from rank i process
                        pi_diskProtonCount,         &       ! how many elements ?
                        pi_mpiDiskProtonType,       &       ! type of elements sent
                        pi_diskProtons,             &       ! stored into array on master only
                        pi_diskProtonCountProcs,    &       ! number of elements stored (master)
                        pi_diskProtonCountOffsets,  &       ! where elements are stored (master)
                        pi_mpiDiskProtonType,       &       ! type of elements stored (master)
                        MASTER_PE,                  &       ! rank of master
                        pi_globalComm,              &       ! communicator handle
                        error                       )       ! error handle
  end if
!
!
!     ...Write out all disk protons on the master processor to disk.
!        Increment new disk proton batch counter on new disk proton file (master only).
!
!
  if (pi_globalMe == MASTER_PE) then

      inquire (unit = pi_diskProtonNewFileID, exist = fileExists)

      if (.not.fileExists) then
           call Driver_abortFlash ("pi_flushDiskProtons2Disk: New disk proton file does not exist!")
      end if

      if (pi_globalDiskProtonCount > 0) then

          write (pi_diskProtonNewFileID)  pi_globalDiskProtonCount
          write (pi_diskProtonNewFileID) (pi_diskProtons (1:PROTON_ATTRCOUNT,1:pi_globalDiskProtonCount))

          pi_diskProtonNewFileNbatches = pi_diskProtonNewFileNbatches + 1

          write (pi_monitorFileUnit,'(a,i8,a,i8,a)') ' flushed ', pi_globalDiskProtonCount, &
               ' Domain Disk   Protons -> New Disk (total # proton batches = ',pi_diskProtonNewFileNbatches,')'
      end if
  end if
!
!
!     ...Ready! 
!
!
  return
end subroutine pi_flushDiskProtons2Disk
