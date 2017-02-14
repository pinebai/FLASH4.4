!!****if* source/flashUtilities/Pipeline/localAPI/pl_localPostRecvMsg
!!
!! NAME
!!  
!!  pl_localPostRecvMsg
!!
!! SYNOPSIS
!! 
!!  call pl_localPostRecvMsg (integer, intent (in) :: channel)
!!
!! DESCRIPTION
!!
!!  Posts a receive message for the specified channel on the local processor.
!!
!! ARGUMENTS
!!
!!  channel : channel index for which to post the receive
!!
!!***

subroutine pl_localPostRecvMsg (channel)

  implicit none

  integer, intent (in) :: channel

  return
end subroutine pl_localPostRecvMsg
