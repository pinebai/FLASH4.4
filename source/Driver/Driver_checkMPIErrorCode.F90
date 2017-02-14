!!****f* source/Driver/Driver_checkMPIErrorCode
!!
!! NAME
!!
!!   Driver_checkMPIErrorCode
!!
!! SYNOPSIS
!!
!!   Driver_checkMPIErrorCode( integer(IN) :: errorCode)
!!
!! ARGUMENTS
!!
!!   errorCode : Return value of an MPI call.
!!
!! DESCRIPTION
!!
!!  Checks the return value of an MPI call.  If all is healthy then
!!  the value is MPI_SUCCESS and the subroutine simply returns.
!!  If we encounter non-MPI_SUCCESS we attempt to learn as much about 
!!  the error as possible.  This includes translating error codes and 
!!  error classes to strings describing the error.  After discovering 
!!  an error we eventually call Driver_abortFlash.
!!
!! NOTES
!!   
!!  When there is an error we print the error message to the main FLASH
!!  logfile.  We are able to pass any size string to Logfile_stampMessage 
!!  because the string argument in Logfile_stampMessage is character (len=*).
!!
!!  "The error codes returned, with the exception of MPI_SUCCESS, are 
!!  defined by each implementation.  This approach allows an MPI
!!  implementation to encode additional data into the error code.  MPI 
!!  also specifies a small set of error classes: integers that divide 
!!  the error codes into a small number of categories." [Ref:Using MPI].
!!
!!  "Several factors limit the ability of MPI calls to return with
!!  meaningful error codes when an error occurs.  MPI may not be able to
!!  detect some errors; other errors may be too expensive to detect in
!!  normal execution mode; finally some errors may be "catastrophic" and
!!  may prevent MPI from returning control to the caller in a consistent
!!  state." [Ref:MPI standard].
!!
!!***

Subroutine Driver_checkMPIErrorCode( errorCode)

implicit none
integer, intent(IN) ::  errorCode

End Subroutine Driver_checkMPIErrorCode
