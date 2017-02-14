!!****if* source/monitors/Logfile/LogfileMain/Logfile_stampMessage
!!
!! NAME
!!
!!  Logfile_stampMessage
!!
!! SYNOPSIS
!!
!!  Logfile_stampMessage(character(len=*)(in)  :: string,
!!                       logical(in),optional :: force )
!!
!! DESCRIPTION
!!   Stamps a simple string to the logfile.
!!
!! ARGUMENTS
!!
!!
!!   string : string to stamp to the logfile
!!
!!   force : if .false., write to logfile only if my PE == MASTER_PE;
!!   force : if .true., logfile is stamped no matter what the PE;
!!           a local logfile (specific to the calling PE) is used,
!!           and the string is prefixed with a timestamp, if PE .ne. MASTER_PE.
!!
!!  NOTES
!!
!!   In general, force is .true. only when called from Driver_abortFlash
!!   Setting force=true can cause very bad performance in multi-processor runs.
!!
!!   MASTER_PE is defined in constants.h .
!!
!!***

#include "Flash.h"

#ifndef FLASH_FLUSH
! Do not actually flush if FLASH_FLUSH is not already defined
#define FLASH_FLUSH(u)
#endif

subroutine Logfile_stampMessage( string,force)

  use Logfile_data, ONLY : log_fileOpen, log_fileOpenLocal, &
       log_keepOpenAfterStamp, log_flushLevel,              &
       log_lun, log_lunLocal, log_globalMe
  use Logfile_interface, ONLY : Logfile_close, Logfile_open, &
       Logfile_getDateTimeStr

  implicit none

#include "constants.h"
  character(len=*), intent(in)           :: string
  logical, intent(in), optional          :: force

  character(len=28)                      :: dateTimeStr
  logical                                :: forceLocal 
  logical                                :: fileIsOpen
  integer :: logUnit
  logical :: logUnitLocal

  forceLocal = .false.
  if (present(force)) forceLocal = force

  ! only master processor writes to logfile
  if ((.not. forceLocal) .and. (log_globalMe .ne. MASTER_PE))  return

  if (forceLocal .and. (log_globalMe .ne. MASTER_PE)) then
     logUnitLocal = .TRUE.
     fileIsOpen   = log_fileOpenLocal
     logUnit      = log_lunLocal
  else
     logUnitLocal = .FALSE.
     fileIsOpen   = log_fileOpen
     logUnit      = log_lun
  end if


  if (.not. fileIsOpen) then
     call Logfile_open(logUnit,logUnitLocal)
     if (forceLocal .and. (log_globalMe .ne. MASTER_PE)) then
        fileIsOpen = log_fileOpenLocal
     else
        fileIsOpen = log_fileOpen
     end if
  end if
  
  if (fileIsOpen) then
     if (logUnitLocal) then
        call Logfile_getDateTimeStr(dateTimeStr)
100     format(" ",A," ",A)
        write (logUnit, 100) dateTimeStr, trim(string)
     else
        write (logUnit, *) trim(string)
     end if
  end if


  if (fileIsOpen .AND.  (logUnitLocal .OR. .NOT. log_keepOpenAfterStamp)) then
     call Logfile_close(logUnitLocal)
  else if (fileIsOpen .AND. log_flushLevel > 0) then
     FLASH_FLUSH(logUnit)
  end if

  return
  
end subroutine Logfile_stampMessage

