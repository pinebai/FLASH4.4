!!****if* source/monitors/Logfile/LogfileMain/Logfile_close
!! NAME
!!   Logfile_close
!!
!! SYNOPSIS
!!   Logfile_close(logical(IN) :: logUnitLocal) 
!!
!! DESCRIPTION
!!   Close the log file
!!
!! ARGUMENTS
!!   logUnitLocal - indicates whether to close the local or global logfile
!!
!! NOTES
!!  variables that begin with "log_" are defined in the fortran 
!!  module Logfile_data.  The prefix "log_" is meant to indicate
!!  that these variables have Logfile unit scope.  Other variables
!!  are local to the individual subroutines
!!
!!***

subroutine Logfile_close(logUnitLocal)

  use Logfile_data, ONLY : log_globalMe, log_fileOpen, log_lun,log_fileOpenLocal, log_lunLocal

#include "constants.h"
  
  implicit none

  logical, optional, intent(IN) :: logUnitLocal


  ! Close "local" (PE-specific) logfile instead of the global one if so requested
  if(present(logUnitLocal)) then
     if(logUnitLocal) then
        close(log_lunLocal)
        log_fileOpenLocal=.false.
        return
     end if
  end if
  
  ! We get here only if logUnitLocal was not present or if it was .FALSE.
  if(.not. log_fileOpen) return
   
  if(log_globalMe == MASTER_PE) close(log_lun)

  log_fileOpen = .false.

end subroutine Logfile_close

