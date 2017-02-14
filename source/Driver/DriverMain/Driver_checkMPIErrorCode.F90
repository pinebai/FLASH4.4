!!****if* source/Driver/DriverMain/Driver_checkMPIErrorCode
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

  use Driver_interface, ONLY : Driver_abortFlash
  use Driver_data, ONLY : dr_globalMe
  implicit none
  include "Flash_mpi.h"
  
  integer, intent(IN) :: errorCode
  
  !Converting integers to strings.
  integer, parameter :: NumConvChars = 20
  character (len=*), parameter :: FMTstrInt = "(i20)"  !Equal to NumConvChars.
  character (len=NumConvChars) :: sIntPE, sInt  !Strings for holding Ints.
  
  !Creating a prefix string consisting of processor identifier.
  character (len=*), parameter :: strPE = "Processor:"
  integer, parameter :: PEDescription = len(strPE) + NumConvChars
  
  !Error messages.
  character (len=MPI_MAX_ERROR_STRING) :: errorString
  character (len=*), parameter :: msgPart1 = &
       " MPI error code outside valid range.  MPI error code:"
  character (len=*), parameter :: msgPart2 = &
       " MPI error (code / class):"
  character (len=*), parameter :: msgPart3 = &
       " MPI error (code / class) corresponds to string:"
  character (len=*), parameter :: msgPart4 = &
       " MPI query call failed (Ignore!).  MPI query error code:"
  
  integer, parameter :: msgLen1 = len(msgPart1) + NumConvChars
  integer, parameter :: msgLen2 = len(msgPart2) + NumConvChars
  integer, parameter :: msgLen3 = len(msgPart3) + MPI_MAX_ERROR_STRING
  integer, parameter :: msgLen4 = len(msgPart4) + NumConvChars
  integer, parameter :: msgLen = PEDescription + &
       max(msgLen1, msgLen2, msgLen3, msgLen4)
  character (len=msgLen) :: msg
  
  !Misc.
  integer :: ierr, stringLen, errorClass
  logical, parameter :: forceWrite = .true.
  
  
  if (errorCode == MPI_SUCCESS) then
     return  !Yeah.  No error.
  end if
  
  
  write(sIntPE,FMTstrInt) dr_globalMe
  write(sInt,FMTstrInt) errorCode
  
  if ((errorCode < MPI_SUCCESS) .or. (errorCode > MPI_ERR_LASTCODE)) then
     msg = strPE // trim(sIntPE) // msgPart1 // trim(sInt)
     print *, msg
     call Driver_abortFlash("[Driver_checkMPIErrorCode]: Invalid MPI error code.")
  else !Print error code.
     msg = strPE // trim(sIntPE) // msgPart2 // trim(sInt)
     print *, msg
  end if
  
  
  
  !Get the description string corresponding to the orginal error code.
  call MPI_Error_string(errorCode, errorString, stringLen, ierr)
  if (ierr == MPI_SUCCESS) then  !Print string.
     msg = strPE // trim(sIntPE) // msgPart3 // errorString(1:stringLen)
  else !Print query error code (probably not useful).
     write(sInt,FMTstrInt) ierr
     msg = strPE // trim(sIntPE) // msgPart4 // trim(sInt)
  end if
  print *, msg
  
  
  
  !Find the more general error code given by the class.
  call MPI_Error_class(errorCode, errorClass, ierr)
  if (ierr == MPI_SUCCESS) then  !Print error class code.
     write(sInt,FMTstrInt) errorClass
     msg = strPE // trim(sIntPE) // msgPart2 // trim(sInt)   
  else !Print query error code (probably not useful).
     write(sInt,FMTstrInt) ierr
     msg = strPE // trim(sIntPE) // msgPart4 // trim(sInt)
  end if
  print *, msg
  
  
  
  !Get the description string corresponding to the class error code.
  call MPI_Error_string(errorClass, errorString, stringLen, ierr)
  if (ierr == MPI_SUCCESS) then  !Print string.
     msg = strPE // trim(sIntPE) // msgPart3 // errorString(1:stringLen)
  else !Print query error code (probably not useful).
     write(sInt,FMTstrInt) ierr
     msg = strPE // trim(sIntPE) // msgPart4 // trim(sInt)
  end if
  print *, msg
  
  
  !Don't attempt to recover from error, just crash.
  call Driver_abortFlash("[Driver_checkMPIErrorCode]: MPI error encountered.")
  
End Subroutine Driver_checkMPIErrorCode
