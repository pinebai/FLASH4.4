!!****if* source/flashUtilities/Pipeline/PipelineMain/Pipeline_localSendItem
!!
!! NAME
!!  
!!  Pipeline_localSendItem
!!
!! SYNOPSIS
!! 
!!  call Pipeline_localSendItem (real,    intent (in)  :: userItem (:),
!!                               integer, intent (in)  :: userProcID,
!!                               logical, intent (out) :: isHandled)
!!
!! DESCRIPTION
!!
!!  The routine will try to add the user supplied item to the local send buffer on
!!  the pipeline channel corresponding to the user supplied channel processor procID.
!!  This can be done, if there is no pending send on the send buffer and if there
!!  is enough space on the send buffer. This part of the routine can be considered
!!  the local feeder of the pipeline.
!!
!! ARGUMENTS
!!
!!  userItem   : the item (array of elements) to be added to the send buffer
!!  userProcID : channel processor ID
!!  isHandled  : is true, if the item was successfully added to the buffer
!!
!! NOTES
!!
!!  If the item supplied by the user contains more elements than the pipeline items
!!  can handle, the program aborts.
!!
!!  It is the user's responsibility to make sure that the supplied channel processor
!!  ID actually belongs to a channel of the pipeline on the local processor. If the
!!  userProcID channel is not found within all the channels of the pipeline on the
!!  local processor, then currently we signal that something must have gone wrong.
!!  Either the pipeline for the local processor was set up incorrectly or userProcID
!!  does not correspond to a channel processor. In this case the only solution is to
!!  abort the program.
!!
!!  In the future we might just simply ignore the sending item request, if userProcID
!!  is not one of the channels. We simply return isHandled = false and continue with
!!  the program. 
!!
!!  Currently the routine is written in such a way that the item will always be added
!!  to the send buffer, i.e. isHandled should always be true. In case there is a
!!  pending send of the buffer, we call the progress sending communication repeatedly
!!  until the send for the particular channel has completed (i.e. the send buffer is
!!  free to be reused). A return of isHandled equal to false indicates an abnormality
!!  and should be catched by the calling routine.
!!
!!***

subroutine Pipeline_localSendItem (userItem, userProcID,   isHandled)

  use Pipeline_data,     ONLY : pl_channelSize, &
                                pl_doLog,       &
                                pl_itemSize,    &
                                pl_logUnit,     &
                                pl_numChannels, &
                                pl_procList,    &
                                pl_sendBuf,     &
                                pl_sendCount,   &
                                pl_sendRequest, &
                                pl_sendStatus


  use Driver_interface,  ONLY : Driver_abortFlash

  use pl_interface,      ONLY : pl_localPostSendMsg

  implicit none

  include "Flash_mpi.h"

  real,    intent (in)  :: userItem (:)
  integer, intent (in)  :: userProcID
  logical, intent (out) :: isHandled

  integer :: channel, itemsCount
  integer :: error
  integer :: i

  integer, parameter :: notFound = -1
!
!
!     ...Check, if the pipeline can handle the user supplied item.
!
!
  if (size (userItem) > pl_itemSize) then
      call Driver_abortFlash ('[Pipeline_localSendItem] ERROR: Size of user supplied item to large!')
  end if
!
!
!     ...Identify the channel on the pipeline corresponding to the user supplied
!        processor ID. Abort run if not found.
!
!
  channel = notFound

  do i = 1, pl_numChannels                        ! It may be necessary to change the
     if (pl_procList (i) == userProcID) then      ! pl_procList data structure to make
         channel = i                              ! the lookup faster (maybe pre-sorting
         exit                                     ! the list)
     end if
  end do

  if (channel == notFound) then
      call Driver_abortFlash ('[Pipeline_localSendItem] ERROR: channel for item not found!')
  end if
!
!
!     ...If there is a pending send in our desired channel, we have to wait until that
!        send completes.
!
!
  if (pl_sendRequest (channel) /= MPI_REQUEST_NULL) then

      call MPI_Wait (pl_sendRequest  (channel), &
                     pl_sendStatus (:,channel), &
                     error                      )

      pl_sendCount (channel) = 0      ! reset the channel item count to 0

      if (pl_doLog) then
          write (pl_logUnit,'(a,i6)') ' Completed send           to proc ID ',userProcID
      end if

  end if
!
!
!     ...We can safely add items to the send buffer if there is no pending send.
!        If the channel is full, post a sending message from that channel.
!
!
  itemsCount = pl_sendCount (channel) + 1

  if (itemsCount > pl_channelSize) then
      call Driver_abortFlash ('[Pipeline_localSendItem] ERROR: channel overflow!')
  end if

  pl_sendBuf   (:, itemsCount , channel) = userItem (:)
  pl_sendCount (                channel) = itemsCount
          
  if (itemsCount == pl_channelSize) then
      call pl_localPostSendMsg (channel)
  end if

  isHandled = .true.
!
!
!    ...Ready!
!
!
  return
end subroutine Pipeline_localSendItem
