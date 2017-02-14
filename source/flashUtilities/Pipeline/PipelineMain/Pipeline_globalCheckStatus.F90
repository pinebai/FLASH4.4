!!****if* source/flashUtilities/Pipeline/PipelineMain/Pipeline_globalCheckStatus
!!
!! NAME
!!  
!!  Pipeline_globalCheckStatus
!!
!! SYNOPSIS
!! 
!!  call Pipeline_globalCheckStatus (logical, intent (out) :: finishedComm,
!!                                   logical, intent (out) :: emptyRecvBuffers,
!!                                   logical, intent (out) :: emptyItemBuffers)
!!
!! DESCRIPTION
!!
!!  Returns the global status of the pipeline. This is a blocking operation that
!!  currently evaluates the following statuses: 1) has all communication inside
!!  the pipeline been performed?, 2) are all receive buffers on all channels on
!!  all processors empty? and 3) are all item buffers on all channels empty?.
!!
!!  Procedure:
!!
!!  The status vector on all processors are all summation reduced and redistributed
!!  on all processors. On each processor the communication entry contains a positive
!!  integer +n for all n sends posted and a collection of negative integers {-m} for
!!  all channel processor ID's on which -m receives completed. Reduction under summation
!!  will give the overall active sends that have not completed. All entries in the
!!  reduced status vector must hence be >= 0. A negative entry indicates a bug someplace.
!!
!!  The receive and item buffer information in the status vector is non-overlapping
!!  and specific to each processor (with entries of +1 or 0). Hence reduction under
!!  summation and redistribution will not change the info contained in the status vector.
!!
!! ARGUMENTS
!!
!!  finishedComm     : is set true, if no more communication is scheduled globally
!!  emptyRecvBuffers : is set true, if all receive buffers on all processors are empty
!!  emptyItemBuffers : is set true, if all item buffers on all processors are empty
!!
!! NOTES
!!
!!  This is a collective operation, hence all processors must call it.
!!
!!***

subroutine Pipeline_globalCheckStatus (finishedComm, emptyRecvBuffers, emptyItemBuffers)

  use Pipeline_data, ONLY : pl_comm,             &
                            pl_doLog,            &
                            pl_logUnit,          &
                            pl_procStatusGlobal, &
                            pl_procStatusLocal,  &
                            pl_size

  implicit none

#include "Pipeline.h"
 include "Flash_mpi.h"

  logical, intent (out) :: finishedComm
  logical, intent (out) :: emptyRecvBuffers
  logical, intent (out) :: emptyItemBuffers

  integer :: error
!
!
!     ...Perform the local status reduction/distribution and asssemble the global
!        status vector. Reset the local status vector to local.
!
!
  call MPI_Allreduce (MPI_IN_PLACE,             &
                      pl_procStatusLocal,       &
                      pl_size * PL_STATUS_SIZE, &
                      FLASH_INTEGER,            &
                      MPI_Sum,                  &
                      pl_comm,                  &
                      error                     )

  pl_procStatusGlobal (:,PL_STATUS_COMM) =   pl_procStatusGlobal (:,PL_STATUS_COMM) &
                                           + pl_procStatusLocal  (:,PL_STATUS_COMM)
  pl_procStatusGlobal (:,PL_STATUS_RECV) =   pl_procStatusLocal  (:,PL_STATUS_RECV)
  pl_procStatusGlobal (:,PL_STATUS_ITEM) =   pl_procStatusLocal  (:,PL_STATUS_ITEM)

  pl_procStatusLocal  (:,PL_STATUS_COMM) = 0
  pl_procStatusLocal  (:,PL_STATUS_RECV) = 0
  pl_procStatusLocal  (:,PL_STATUS_ITEM) = 0

  if (pl_doLog) then
      call pl_printGlobalStatusVector (pl_logUnit)
  end if

  finishedComm     = all (pl_procStatusGlobal (:,PL_STATUS_COMM) == 0)
  emptyRecvBuffers = all (pl_procStatusGlobal (:,PL_STATUS_RECV) == 0)
  emptyItemBuffers = all (pl_procStatusGlobal (:,PL_STATUS_ITEM) == 0)

  if (pl_doLog) then
      write (pl_logUnit,*) ' finishedComm     = ',finishedComm
      write (pl_logUnit,*) ' emptyRecvBuffers = ',emptyRecvBuffers
      write (pl_logUnit,*) ' emptyItemBuffers = ',emptyItemBuffers
  end if
!
!
!    ...Ready!
!
!
  return
end subroutine Pipeline_globalCheckStatus
