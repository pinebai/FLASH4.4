!!****f* source/Driver/Driver_abortFlash
!!
!! NAME
!!
!!  Driver_abortFlash
!!
!! SYNOPSIS
!!
!!  Driver_abortFlash(character(len=*)(IN) :: errorMessage)
!!
!! DESCRIPTION
!!
!!  Write an error message to the logfile and abort FLASH.
!!  Attempts to shut down all processes (using MPI_Abort()).
!!  If you wish to call Driver_abortFlash from a 'c' routine
!!  use the API routine Driver_abortFlashC
!!
!! ARGUMENTS
!!
!!  errorMessage :    A string to write to the logfile (presumably 
!!                    indicating what went wrong).
!!
!! NOTES
!!
!!  This function's implementation never returns control to the caller.
!!  
!!
!!***

subroutine Driver_abortFlash (errorMessage)
  
  implicit none

  character(len=*), intent(in) :: errorMessage

  return
end subroutine Driver_abortFlash
