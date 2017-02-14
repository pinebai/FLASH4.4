!!****f* source/Driver/Driver_logMemoryUsage
!!
!! NAME
!!
!!  Driver_logMemoryUsage
!!
!! SYNOPSIS
!!
!!  Driver_logMemoryUsage(character(len=*)(IN) :: callsite)
!!
!! DESCRIPTION
!!
!!  Logs memory usage
!!
!!
!! ARGUMENTS
!!
!!  callsite -    A string storing from where we call this subroutine.
!!
!! NOTES
!!
!!  This routine is a no-op if the preprocessor symbol FLASH_USE_MEMORYUSAGE is undefined.
!!***

subroutine Driver_logMemoryUsage (callsite)
  implicit none
  character(len=*), intent(IN) :: callsite
end subroutine Driver_logMemoryUsage
