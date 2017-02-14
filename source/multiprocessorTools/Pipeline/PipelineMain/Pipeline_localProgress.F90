!!****if* source/multiprocessorTools/Pipeline/PipelineMain/Pipeline_localProgress
!!
!! NAME
!!  
!!  Pipeline_localProgress
!!
!! SYNOPSIS
!! 
!!  call Pipeline_localProgress ()
!!
!! DESCRIPTION
!!
!!  This routine can be considered the motor of the pipeline. It makes sure the items
!!  are transported through the channels on the local processor. As such the routine
!!  should be called often to ensure best possible communication progress between the
!!  channels. Progress on a processor can stall, if the item buffer is full.
!!
!! ARGUMENTS
!!
!!  none
!!
!! NOTES
!!
!!  All MPI ranks must call this subroutine to avoid possible deadlock. This subroutine
!!  must remain non-blocking.
!!
!!***

subroutine Pipeline_localProgress ()

  use Pipeline_data,     ONLY : pl_doLog,              &
                                pl_itemSize,           &
                                pl_logUnit,            &
                                pl_numChannels,        &
                                pl_procList,           &
                                pl_rank,               &
                                pl_recvCount,          &
                                pl_recvIndex,          &
                                pl_recvRequest,        &
                                pl_recvStatus,         &
                                pl_sendCount,          &
                                pl_sendIndex,          &
                                pl_sendRequest,        &
                                pl_sendStatus

  use Driver_interface,  ONLY : Driver_abortFlash,       &
                                Driver_checkMPIErrorCode

  use pl_interface,      ONLY : pl_localPostRecvMsg,  &
                                pl_localSaveRecvItems

  implicit none

#include "Pipeline.h"
 include "Flash_mpi.h"

  logical :: isSaved

  integer :: n, nCompleted
  integer :: error
  integer :: channel
  integer :: nItems, nReals
  integer :: procID
!
!
!     ...Determines progress of receive communications inside the pipeline on the
!        local processor. Receiving channels that have received items, but could not
!        store them into the items buffer due to buffer overflow, are identified
!        and the corresponding receive buffers are emptied (items moved to the items
!        buffer). Next all completed receive channels are identified and their
!        buffers are processed. Once the receive buffers are empty, new speculative
!        receives are posted.
!
!
  if (pl_numChannels > 0) then

      do channel = 1, pl_numChannels

         if (pl_recvRequest (channel) == MPI_REQUEST_NULL .and. &
             pl_recvCount   (channel)  > 0                      ) then

             call pl_localSaveRecvItems   (channel,  isSaved)

             if (isSaved) then
                 call pl_localPostRecvMsg (channel)
             end if
         end if
      end do
!
!
!     ...Test all receive channels for new messages.  Save the corresponding items
!        (if possible) and then post a new receive.
!
!
      call MPI_Testsome (pl_numChannels, &    ! how many should be tested
                         pl_recvRequest, &    ! handles to be tested
                         nCompleted,     &    ! number of completed requests
                         pl_recvIndex,   &    ! 1d array of indices of operations that completed
                         pl_recvStatus,  &    ! 2d array of status objects for operations that completed
                         error           )

      call Driver_checkMPIErrorCode (error)

      do n = 1, nCompleted

         channel = pl_recvIndex  (n)              ! get the channel of the completed receive
         procID  = pl_recvStatus (MPI_SOURCE,n)   ! get the procID from the source handle in status field

         if (procID /= pl_procList (channel)) then
             call Driver_abortFlash ('[Pipeline_localProgress] ERROR: procID mismatch!')
         end if

         call MPI_Get_count (pl_recvStatus (:,n), &   ! use the handles of the status field
                             FLASH_REAL,          &   ! to find out how many reals were
                             nReals,              &   ! received from the completed receive
                             error                )   ! operation

         call Driver_checkMPIErrorCode (error)

         if (mod (nReals, pl_itemSize) > 0) then
             call Driver_abortFlash ('[Pipeline_localProgress] ERROR: # reals / item size mismatch!')
         end if

         nItems = nReals / pl_itemSize                ! each item can contains many elements (reals)
         pl_recvCount (channel) = nItems              ! store number of items received in channel

         if (pl_doLog) then
             write (pl_logUnit,'(2(a,i6))') ' Received  ',nItems,' items from proc ID ',procID
         end if

         if (nItems > 0) then
             call pl_localSaveRecvItems   (channel,  isSaved)
             if (isSaved) then
                 call pl_localPostRecvMsg (channel)
             end if
         end if

      end do
!
!
!     ...Test all sending channels for completed sends. Those that have completed
!        will have their request status reset to MPI_REQUEST_NULL and have their
!        buffer counts reset to zero.
!
!
      call MPI_Testsome (pl_numChannels, &    ! how many should be tested
                         pl_sendRequest, &    ! handles to be tested
                         nCompleted,     &    ! number of completed requests
                         pl_sendIndex,   &    ! 1d array of indices of operations that completed
                         pl_sendStatus,  &    ! 2d array of status objects for operations that completed
                         error           )

      call Driver_checkMPIErrorCode (error)

      do n = 1, nCompleted

         channel = pl_sendIndex  (n)          ! get the channel of the completed send
         procID  = pl_procList (channel)      ! the procID to where the send has been done
         pl_sendCount (channel) = 0           ! reset the channel item count to 0

         if (pl_doLog) then
             write (pl_logUnit,'(a,i6)') ' Completed send           to proc ID ',procID
         end if

      end do

  end if
!
!
!    ...Ready!
!
!
  return
end subroutine Pipeline_localProgress
