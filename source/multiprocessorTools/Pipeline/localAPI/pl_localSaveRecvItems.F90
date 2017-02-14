!!****if* source/multiprocessorTools/Pipeline/localAPI/pl_localSaveRecvItems
!!
!! NAME
!!  
!!  pl_localSaveRecvItems
!!
!! SYNOPSIS
!! 
!!  call pl_localSaveRecvItems (integer, intent (in)  :: channel,
!!                              logical, intent (out) :: isSaved)
!!
!! DESCRIPTION
!!
!!  Saves all received items in the specified channel on the local processor.
!!
!! ARGUMENTS
!!
!!  channel : channel index for which to save the received items
!!  isSaved : if true, all items have been successfully saved
!!
!!***

subroutine pl_localSaveRecvItems (channel, isSaved)

  implicit none

  integer, intent (in)  :: channel
  logical, intent (out) :: isSaved

  isSaved = .false.

  return
end subroutine pl_localSaveRecvItems
