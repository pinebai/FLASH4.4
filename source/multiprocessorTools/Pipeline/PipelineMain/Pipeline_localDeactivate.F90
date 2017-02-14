!!****if* source/multiprocessorTools/Pipeline/PipelineMain/Pipeline_localDeactivate
!!
!! NAME
!!  
!!  Pipeline_localDeactivate
!!
!! SYNOPSIS
!! 
!!  call Pipeline_localDeactivate (logical, optional, intent (in) :: doAsyncReturn)
!!
!! DESCRIPTION
!!
!!  Deactivates the local pipeline environment. Deactivation means smooth shutdown
!!  of the local pipeline communitation environment. There are two ways the local
!!  pipeline communicator can shut down: 1) all processors return at will or 2) all
!!  processors return at the same time.
!!
!! ARGUMENTS
!!
!!  doAsyncReturn : if set true, each processor returns at will
!!
!! NOTES
!!
!!  This routine should be called at the end of a pipeline application. Deactivation
!!  is essential before destroying the pipeline. One should never destruct a pipeline
!!  before deactivating it.
!!
!!***

subroutine Pipeline_localDeactivate (doAsyncReturn)

  use Pipeline_data, ONLY : pl_comm,              &
                            pl_doLog,             &
                            pl_logUnit,           &
                            pl_numChannels,       &
                            pl_pipelineActive,    &
                            pl_rank,              &
                            pl_recvRequest,       &
                            pl_recvStatus,        &
                            pl_sendRequest,       &
                            pl_sendStatus,        &
                            pl_size

  implicit none

  include "Flash_mpi.h"

  logical, optional, intent (in) :: doAsyncReturn

  logical :: doSyncReturn

  integer :: channel
  integer :: error
  integer :: handle
!
!
!     ...If pipeline is currently active and there are channels present, force all
!        local sends and receives to finish.
!
!
  if (pl_pipelineActive) then

      if (pl_numChannels > 0) then

          do channel = 1, pl_numChannels
             handle = pl_recvRequest (channel)
             if (handle /= MPI_REQUEST_NULL) then
                 call MPI_Cancel       (handle, error)
                 call Driver_checkMPIErrorCode (error)
             end if
          end do

          call MPI_Waitall (pl_numChannels, &
                            pl_recvRequest, &
                            pl_recvStatus,  &
                            error           )

          call Driver_checkMPIErrorCode (error)

          do channel = 1, pl_numChannels
             handle = pl_sendRequest (channel)
             if (handle /= MPI_REQUEST_NULL) then
                 call MPI_Cancel       (handle, error)
                 call Driver_checkMPIErrorCode (error)
             end if
          end do

          call MPI_Waitall (pl_numChannels, &
                            pl_sendRequest, &
                            pl_sendStatus,  &
                            error           )

          call Driver_checkMPIErrorCode (error)

      end if

      if (pl_size > 1) then
          if (present (doAsyncReturn)) then
              doSyncReturn = .not. doAsyncReturn
          else
              doSyncReturn = .true.
          end if
          if (doSyncReturn) then
              call MPI_Barrier     (pl_comm, error)
              call Driver_checkMPIErrorCode (error)
          end if
      end if

  end if
!
!
!     ...The pipeline is considered now deactivated.
!
!    
  pl_pipelineActive = .false.

  if (pl_doLog) then
      write (pl_logUnit,'(a,i6)') ' Pipeline deactived       on proc ID ', pl_rank
  end if
!
!
!    ...Ready!
!
!
  return
end subroutine Pipeline_localDeactivate
