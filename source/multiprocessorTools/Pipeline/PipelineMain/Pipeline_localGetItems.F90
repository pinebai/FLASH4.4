!!****if* source/multiprocessorTools/Pipeline/PipelineMain/Pipeline_localGetItems
!!
!! NAME
!!  
!!  Pipeline_localGetItems
!!
!! SYNOPSIS
!! 
!!  call Pipeline_localGetItems (real,    intent (inout) :: userArray (:,:),
!!                               integer, intent (out)   :: userCount)
!!
!! DESCRIPTION
!!
!!  The routine copies items currently sitting in the item buffer of the local
!!  processor to the user supplied array. Once the items in the buffer have been
!!  copied, the item buffer counter is reset. This routine constitutes thus the
!!  local retriever of the pipeline.
!!
!! ARGUMENTS
!!
!!  userArray : the user supplied array to hold the copied items
!!  userCount : the number of items copied
!!
!! NOTES
!!
!!  The user supplied array must be able to hold all elements of each item. If this
!!  is not the case, the program aborts.
!!
!!***

subroutine Pipeline_localGetItems (userArray, userCount)

  use Pipeline_data,     ONLY : pl_doLog,            &
                                pl_itemBuf,          &
                                pl_itemCount,        &
                                pl_itemSize,         &
                                pl_logUnit,          &
                                pl_procItemBalance,  &
                                pl_rank

  use Driver_interface,  ONLY : Driver_abortFlash

#include "Pipeline.h"

  implicit none

  real,    intent (inout) :: userArray (:,:)
  integer, intent (out)   :: userCount

  integer :: firstItem
  integer :: itemsToCopy
!
!
!     ...Check, if number of elements per item fits into the user array.
!
!
  if (size (userArray,1) < pl_itemSize) then
      call Driver_abortFlash ('[Pipeline_localGetItems] ERROR: User array cannot hold item elements!')
  end if
!
!
!     ...Decide how many items to copy. Always copy from the end of the items buffer, to
!        avoid necessary reshuffling of remaining items in buffer.
!
!
  itemsToCopy = min (pl_itemCount, size (userArray,2))

  if (itemsToCopy > 0) then
      firstItem = pl_itemCount - itemsToCopy + 1
      userArray (:,1:itemsToCopy) = pl_itemBuf (:,firstItem:pl_itemCount)
      pl_itemCount = pl_itemCount - itemsToCopy
  end if

  userCount = itemsToCopy

  if (pl_doLog .and. userCount > 0) then
      write (pl_logUnit,'(a,i6,a)') ' Copied                              ',userCount, &
                                    ' buffer items to user supplied array '
  end if
!
!
!     ...Update item balance on current processor.
!
!
  pl_procItemBalance = pl_procItemBalance - userCount
!
!
!    ...Ready!
!
!
  return
end subroutine Pipeline_localGetItems
