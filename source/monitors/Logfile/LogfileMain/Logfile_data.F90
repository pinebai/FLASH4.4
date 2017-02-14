!!****if* source/monitors/Logfile/LogfileMain/Logfile_data
!!
!! NAME
!!  Logfile_data
!!
!! SYNOPSIS
!!
!!   use Logfile_data
!!
!! DESCRIPTION
!!
!!  Holds the data needed by the Logfile Unit
!!
!!
!! CHANGES
!! 
!!  Removed log_restart         2008-02-08 KW
!!  Logic for keeping logs open April 2018 KW
!!
!!***

module Logfile_data

#include "constants.h"
#include "Flash.h"

  !! Unit number and file name for the global logfile
  integer, parameter                     :: log_lun   = 25 ! Logical unit number for file
  character(len=MAX_STRING_LENGTH), save :: log_fileName  = "/dev/stderr"  ! Name of log file
  logical, save                          :: log_fileOpen = .false.
  logical, save                          :: log_keepOpenAfterStamp = .FALSE.
  integer, save                          :: log_flushLevel = 1 ! maybe try to flush after messages, if this is > 0

  !! Unit number and filename for the logfile local to each processor
  integer, parameter                     :: log_lunLocal = 26
  character(len=MAX_STRING_LENGTH), save :: log_fileNameLocal  ! per proc log file
  logical, save                          :: log_fileOpenLocal=.false.

  !DEV: may want this eventually for restart with only flash.par
  !character(len=MAX_STRING_LENGTH), save :: parmFile
  character(len=MAX_STRING_LENGTH), save :: log_runNum
  character(len=MAX_STRING_LENGTH), save :: log_runComment
  character, parameter                   :: log_endOfLine = ' '

  character(len=MAX_STRING_LENGTH), save, allocatable, dimension(:,:) :: log_strArr

  logical, save :: log_enableGcMaskLogging = .false.
  integer, save :: log_globalMe,log_globalNumProcs, log_globalComm

end module Logfile_data
