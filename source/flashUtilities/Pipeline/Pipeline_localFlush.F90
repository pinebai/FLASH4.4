!!****f* source/flashUtilities/Pipeline/Pipeline_localFlush
!!
!! NAME
!!  
!!  Pipeline_localFlush
!!
!! SYNOPSIS
!! 
!!  call Pipeline_localFlush (logical, intent (in) :: fullestChannelOnly)
!!
!! DESCRIPTION
!!
!!  Flushes the local pipeline processor. When calling this routine, either
!!  only the fullest channel will be sent (fullestChannelOnly = .true.) or
!!  all channels with non-empty sending buffers are forced to send.
!!
!! ARGUMENTS
!!
!!  fullestChannelOnly : if true, only the fullest channel will send
!!
!! NOTES
!!
!!  All processors need to call this routine in order to move items through the
!!  pipeline, especially near the end of the application.
!!
!!***

subroutine Pipeline_localFlush (fullestChannelOnly)

  implicit none

  logical, intent (in) :: fullestChannelOnly

  return
end subroutine Pipeline_localFlush
