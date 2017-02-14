!!****if* source/multiprocessorTools/Pipeline/PipelineMain/pl_localSaveRecvItems
!!
!! NAME
!!  
!!  pl_localSaveRecvItems
!!
!! SYNOPSIS
!! 
!!  call pl_localSaveRecvItems (integer, intent (in)  :: channel,
!!                             logical, intent (out) :: isSaved)
!!
!! DESCRIPTION
!!
!!  Saves received items in the specified channel on the local processor.
!!  Three situations can arise here:
!!
!!    1) all received items in the receive buffer can be saved to
!!       the items buffer.
!!
!!    2) only a fraction of the received items in the receive buffer
!!       can be saved to the items buffer, because of items buffer
!!       overflow.
!!
!!    3) no received items can be saved, because the items buffer
!!       is currently full.
!!
!!  Situations 2) and 3) are common, if the pipeline was created with
!!  channel sizes greater than maximum items buffer sizes. Although
!!  the pipeline will continue to transport items, it will do so not
!!  very efficiently. It is advised to declare maximum items buffer
!!  sizes larger than channel sizes on each processor.
!!
!! ARGUMENTS
!!
!!  channel : channel index for which to save the received items
!!  isSaved : if true, all items have been successfully saved
!!
!!***

subroutine pl_localSaveRecvItems (channel, isSaved)

  use Pipeline_data,     ONLY : pl_doLog,            &
                                pl_itemBuf,          &
                                pl_itemCount,        &
                                pl_itemSize,         &
                                pl_logUnit,          &
                                pl_maxItems,         &
                                pl_procList,         &
                                pl_rank,             &
                                pl_recvBuf,          &
                                pl_recvCount

#include "Pipeline.h"

  implicit none

  integer, intent (in)  :: channel
  logical, intent (out) :: isSaved

  integer :: bufferBeg, bufferEnd
  integer :: n, nItems, nSave
  integer :: procID
  integer :: recvBeg
!
!
!     ...Save the received items (if any) and write this action to the log file (if requested).
!
!
  nItems = pl_recvCount (channel)

  if (nItems > 0) then

      n = pl_itemSize

      bufferBeg = pl_itemCount + 1
      bufferEnd = pl_itemCount + nItems

      if (bufferEnd <= pl_maxItems) then

          pl_itemBuf (1:n, bufferBeg:bufferEnd) = pl_recvBuf (1:n, 1:nItems, channel)
          pl_itemCount = bufferEnd
          pl_recvCount (channel) = 0

          isSaved = .true.

          if (pl_doLog) then
              procID = pl_procList (channel)
              write (pl_logUnit,'(a,i6,2(a,2(i6)),a)') ' Saved receive items    from proc ID ',procID, & 
                                                       ' : copy from msg slice (',1,nItems,            &
                                                       ' ) to buffer slice (', bufferBeg, bufferEnd, ')'
          end if

      else if (bufferBeg <= pl_maxItems) then

          nSave = pl_maxItems - bufferBeg + 1
          recvBeg = nItems - nSave + 1

          pl_itemBuf (1:n, bufferBeg:pl_maxItems) = pl_recvBuf (1:n, recvBeg:nItems, channel)
          pl_itemCount = pl_maxItems
          pl_recvCount (channel) = nItems - nSave

          isSaved = .false.

          if (pl_doLog) then
              procID = pl_procList (channel)
              write (pl_logUnit,'(a,i6,2(a,2(i6)),a)') ' Saved receive items    from proc ID ',procID, & 
                                                       ' : copy from msg slice (',recvBeg,nItems,      &
                                                       ' ) to buffer slice (', bufferBeg, pl_maxItems, ')'
          end if
      else
          isSaved = .false.
      end if

  end if
!
!
!    ...Ready!
!
!
  return
end subroutine pl_localSaveRecvItems
