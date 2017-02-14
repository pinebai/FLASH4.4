!!****if* source/monitors/Logfile/LogfileMain/Logfile_break
!!
!! NAME
!!  Logfile_break
!!
!! SYNOPSIS
!!
!!  Logfile_break(character(in):: char)
!!
!! DESCRIPTION
!!  Writes a horizontal line of a given char to the Logfile.  This is used to
!!  break up different sections of the Logfile, performance results from runtime stamping
!!  etc.
!! 
!!  ===========================================================================
!!
!! ARGUMENTS
!!
!!  char - single character used to create line in Logfile
!!
!! NOTES
!!  variables that begin with "log_" are defined in the fortran 
!!  module Logfile_data.  The prefix "log_" is meant to indicate
!!  that these variables have Logfile unit scope.  Other variables
!!  are local to the individual subroutines
!! 
!! 
!!
!!***

subroutine Logfile_break ( char)

  use Logfile_data, ONLY : log_fileOpen, log_lun, log_globalMe
  use Logfile_interface, ONLY : Logfile_open

  implicit none

#include "constants.h"


  
  character(len=1), intent(in)                   :: char
  character(len=MAX_STRING_LENGTH-2) :: buff
  
  character(len=1)                   :: lchar
  !logical                            :: do_open = .false.
  integer :: logUnit
  logical :: logUnitLocal=.false.

  if (log_globalMe .ne. MASTER_PE) return  ! only master processor writes to logfile
  
  if (ichar(char) < 32) then
     lchar = ' '
  else
     lchar = char
  end if
  buff = repeat(lchar, MAX_STRING_LENGTH-2) ! 78 characters

  if (.not. log_FileOpen) then
     call Logfile_open(logUnit,logUnitLocal)
  end if


  if (log_fileOpen) then
     write(log_lun, *) buff
  else
     print *, "ERROR: can not write to logfile"
  end if

  return
end subroutine Logfile_break
