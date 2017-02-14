!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/ProtonDetection/pi_flushScreenProtons2Disk
!!
!! NAME
!!
!!  pi_flushScreenProtons2Disk
!!
!! SYNOPSIS
!!
!!  call pi_flushScreenProtons2Disk ()
!!
!! DESCRIPTION
!!
!!  When calling this routine, the screen protons accumulated by each processor are all
!!  send to the master processor and processed for writing out to disk. The following
!!  steps are performed:
!!
!!     1) Gather at the master processor the info of how many screen protons each
!!        processor currently has. Calculate the offsets to prepare for storage
!!        of the screen protons on the master processor.
!!
!!     2) Gather all the screen protons on the master processor. The screen protons
!!        on the master remain in position, the others are appended according to the
!!        offsets calculated previously.
!!
!!     3) On the master processor, sort the collection of all screen protons into
!!        buckets according to their detector number. There are as many buckets
!!        as are screen detectors. If a bucket is full, trigger writing to disk to
!!        the assigned detector printout file.
!!
!!     4) Once all screen protons have been processed and all buckets have been emptied,
!!        reset the screen proton counter to 0 on all processors and reset the bucket
!!        counters on the master processor to 0 as well.
!!
!! ARGUMENTS
!!
!! NOTES
!!
!!  The routine has been designed in such a way that it can be called multiple times
!!  during a proton imaging simulation. The screen protons written out will be simply
!!  appended to the existing printout detector files.
!!
!!***

subroutine pi_flushScreenProtons2Disk ()

  use ProtonImaging_data,  ONLY : pi_detectorFilesID,          &
                                  pi_detectorLNwriteFormat,    &
                                  pi_globalComm,               &
                                  pi_globalMe,                 &
                                  pi_globalNumProcs,           &
                                  pi_globalScreenProtonCount,  &
                                  pi_maxProtonCount,           &
                                  pi_monitorFileUnit,          &
                                  pi_mpiScreenProtonType,      &
                                  pi_numberOfDetectors,        &
                                  pi_screenProtonBuckets,      &
                                  pi_screenProtonBucketCount,  &
                                  pi_screenProtonBucketSize,   &
                                  pi_screenProtonCount,        &
                                  pi_screenProtonCountOffsets, &
                                  pi_screenProtonCountProcs,   &
                                  pi_screenProtonDiagnostics,  &
                                  pi_screenProtonRecordCount,  &
                                  pi_screenProtons
  
  use Driver_interface,    ONLY : Driver_abortFlash

  implicit none

#include "Flash.h"
#include "constants.h"
#include "ProtonImaging.h"

  include "Flash_mpi.h"

  integer :: bucketCount
  integer :: detector
  integer :: error
  integer :: id
  integer :: offset
  integer :: proton
  integer :: rank
  integer :: recordCount

  real    :: JvScreen, KxScreen, KyScreen, KzScreen
  real    :: xScreen, yScreen
!
!
!     ...Gather the number of screen protons from all processors at the master processor.
!
!
  call MPI_Gather (pi_screenProtonCount,         &       ! this is being sent from rank i process
                   1,                            &       ! how many elements ?
                   MPI_INTEGER,                  &       ! type of elements sent
                   pi_screenProtonCountProcs,    &       ! stored into i-th array element on master only
                   1,                            &       ! number of elements stored (master)
                   MPI_INTEGER,                  &       ! type of elements stored (master)
                   MASTER_PE,                    &       ! rank of master
                   pi_globalComm,                &       ! communicator handle
                   error                         )       ! error handle
!
!
!     ...Calculate the offsets on the master only. The offset for the master is set
!        as zero, as we do not want to move screen protons on the master. Check, if
!        the total number of expected screen protons does not exceed the maximum
!        declared storage.
!
!        Next gather all screen protons from all processors on the master processor.
!        The master at this point has already the offset info for each processor.
!        Do not move screen protons on the master.
!
!
  if (pi_globalMe == MASTER_PE) then

      pi_screenProtonCountOffsets (MASTER_PE) = 0

      offset = pi_screenProtonCountProcs (MASTER_PE)

      do rank = 0,MASTER_PE-1
         pi_screenProtonCountOffsets (rank) = offset
         offset = offset + pi_screenProtonCountProcs (rank)
      end do

      do rank = MASTER_PE+1,pi_globalNumProcs-1
         pi_screenProtonCountOffsets (rank) = offset
         offset = offset + pi_screenProtonCountProcs (rank)
      end do

      pi_globalScreenProtonCount = offset   ! this is the new final (total) screen proton count on the master

      if (pi_globalScreenProtonCount > pi_maxProtonCount) then
          call Driver_abortFlash ('[pi_flushScreenProtons2Disk] ERROR: global # of screen protons > storage size')
      end if

      call MPI_Gatherv (MPI_IN_PLACE,                 &       ! do not move screen protons on master
                        pi_screenProtonCount,         &       ! irrelevant
                        pi_mpiScreenProtonType,       &       ! irrelevant
                        pi_screenProtons,             &       ! stored into array on master only
                        pi_screenProtonCountProcs,    &       ! number of elements stored (master)
                        pi_screenProtonCountOffsets,  &       ! where elements are stored (master)
                        pi_mpiScreenProtonType,       &       ! type of elements stored (master)
                        MASTER_PE,                    &       ! rank of master
                        pi_globalComm,                &       ! communicator handle
                        error                         )       ! error handle

      if (pi_globalScreenProtonCount > 0) then
          write (pi_monitorFileUnit,'(a,i8,a)') ' stored  ',pi_globalScreenProtonCount,' Domain Screen Protons'
      end if

  else

      call MPI_Gatherv (pi_screenProtons,             &       ! this is being sent from rank i process
                        pi_screenProtonCount,         &       ! how many elements ?
                        pi_mpiScreenProtonType,       &       ! type of elements sent
                        pi_screenProtons,             &       ! stored into array on master only
                        pi_screenProtonCountProcs,    &       ! number of elements stored (master)
                        pi_screenProtonCountOffsets,  &       ! where elements are stored (master)
                        pi_mpiScreenProtonType,       &       ! type of elements stored (master)
                        MASTER_PE,                    &       ! rank of master
                        pi_globalComm,                &       ! communicator handle
                        error                         )       ! error handle
  end if
!
!
!     ...Process all screen protons (if any) on the master processor. Distribute (sort) screen protons
!        into their respective buckets and empty them when appropriate.
!
!
  if (pi_globalMe == MASTER_PE) then

      if (pi_globalScreenProtonCount == 0) then
          return
      end if

      pi_screenProtonBucketCount (1:pi_numberOfDetectors) = 0

      if (pi_screenProtonDiagnostics) then

          recordCount = 6

          if (recordCount /= pi_screenProtonRecordCount) then
              call Driver_abortFlash ('[pi_flushScreenProtons2Disk] ERROR: bad size of screen proton buckets!')
          end if

          do proton = 1,pi_globalScreenProtonCount

             xScreen  = pi_screenProtons (SCREEN_POSX,proton)
             yScreen  = pi_screenProtons (SCREEN_POSY,proton)
             detector = pi_screenProtons (SCREEN_DETC,proton)
             JvScreen = pi_screenProtons (SCREEN_DGJV,proton)
             KxScreen = pi_screenProtons (SCREEN_DGKX,proton)
             KyScreen = pi_screenProtons (SCREEN_DGKY,proton)
             KzScreen = pi_screenProtons (SCREEN_DGKZ,proton)

             bucketCount = pi_screenProtonBucketCount (detector)
             bucketCount = bucketCount + 1

             pi_screenProtonBuckets (1,bucketCount,detector) = xScreen
             pi_screenProtonBuckets (2,bucketCount,detector) = yScreen
             pi_screenProtonBuckets (3,bucketCount,detector) = JvScreen
             pi_screenProtonBuckets (4,bucketCount,detector) = KxScreen
             pi_screenProtonBuckets (5,bucketCount,detector) = KyScreen
             pi_screenProtonBuckets (6,bucketCount,detector) = KzScreen

             if (bucketCount == pi_screenProtonBucketSize) then
                 id = pi_detectorFilesID (detector)
                 write (id, pi_detectorLNwriteFormat) (pi_screenProtonBuckets (1:recordCount,1:bucketCount,detector))
                 pi_screenProtonBucketCount (detector) = 0
             else
                 pi_screenProtonBucketCount (detector) = bucketCount
             end if

          end do

      else

          recordCount = 2

          if (recordCount /= pi_screenProtonRecordCount) then
              call Driver_abortFlash ('[pi_flushScreenProtons2Disk] ERROR: bad size of screen proton buckets!')
          end if

          do proton = 1,pi_globalScreenProtonCount

             xScreen  = pi_screenProtons (SCREEN_POSX,proton)
             yScreen  = pi_screenProtons (SCREEN_POSY,proton)
             detector = pi_screenProtons (SCREEN_DETC,proton)

             bucketCount = pi_screenProtonBucketCount (detector)
             bucketCount = bucketCount + 1

             pi_screenProtonBuckets (1,bucketCount,detector) = xScreen
             pi_screenProtonBuckets (2,bucketCount,detector) = yScreen

             if (bucketCount == pi_screenProtonBucketSize) then
                 id = pi_detectorFilesID (detector)
                 write (id, pi_detectorLNwriteFormat) (pi_screenProtonBuckets (1:recordCount,1:bucketCount,detector))
                 pi_screenProtonBucketCount (detector) = 0
             else
                 pi_screenProtonBucketCount (detector) = bucketCount
             end if

          end do

      end if
!
!
!     ...Empty remaining unfull buckets (if any).
!
!
      do detector = 1,pi_numberOfDetectors

         bucketCount = pi_screenProtonBucketCount (detector)

         if (bucketCount > 0) then
             id = pi_detectorFilesID (detector)
             write (id, pi_detectorLNwriteFormat) (pi_screenProtonBuckets (1:recordCount,1:bucketCount,detector))
             pi_screenProtonBucketCount (detector) = 0
         end if

      end do

      write (pi_monitorFileUnit,'(a,i8,a)') ' flushed ',pi_globalScreenProtonCount, &
                                            ' Domain Screen Protons -> Detector Files'

  end if
!
!
!     ...Ready! 
!
!
  return
end subroutine pi_flushScreenProtons2Disk
