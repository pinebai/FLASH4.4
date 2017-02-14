!!****f* source/flashUtilities/Pipeline/Pipeline_localGetItems
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

  implicit none

  real,    intent (inout) :: userArray (:,:)
  integer, intent (out)   :: userCount

  userCount = 0
  userArray (:,:) = 0.0

  return
end subroutine Pipeline_localGetItems
