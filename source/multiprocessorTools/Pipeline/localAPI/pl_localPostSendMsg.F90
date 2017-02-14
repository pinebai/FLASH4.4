!!****if* source/multiprocessorTools/Pipeline/localAPI/pl_localPostSendMsg
!!
!! NAME
!!  
!!  pl_localPostSendMsg
!!
!! SYNOPSIS
!! 
!!  call pl_localPostSendMsg (integer, intent (in) :: channel)
!!
!! DESCRIPTION
!!
!!  Posts a sending message for the specified channel on the local processor.
!!
!! ARGUMENTS
!!
!!  channel : channel index for which to post the send
!!
!!***

subroutine pl_localPostSendMsg (channel)

  implicit none

  integer, intent (in) :: channel

  return
end subroutine pl_localPostSendMsg
