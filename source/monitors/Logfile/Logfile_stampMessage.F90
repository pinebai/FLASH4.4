!!****f* source/monitors/Logfile/Logfile_stampMessage
!!
!! NAME
!!
!!  Logfile_stampMessage
!!
!! SYNOPSIS
!!
!!  Logfile_stampMessage(character(len=*)(in)  :: string,
!!                       logical(in),optional :: force)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!
!!   string : message to stamp to the logfile
!!
!!   force : if force=true, the Logfile is stamped no matter what myPE calls
!!           Otherwise, only the MASTER_PE can stamp the logfile
!!
!!  NOTES
!!
!!    In general, only Driver_abortFlash should set force=true.
!!    Otherwise, extremely slow behaviour or fatal errors may occur in multi-
!!    processor runs
!!
!!***


!stamp only a string message

subroutine Logfile_stampMessage( string,force)

 
  implicit none
  character(len=*), intent(in)           :: string
  logical, intent(in), optional          :: force 

  return
  
end subroutine Logfile_stampMessage

