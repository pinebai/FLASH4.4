!!****f* source/monitors/Logfile/Logfile_break
!!
!! NAME
!!  Logfile_break
!!
!! SYNOPSIS
!!  Logfile_break(character(in) :: char)
!!
!! DESCRIPTION
!!
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
!! 
!! 
!! 
!!
!!***

subroutine Logfile_break ( char)

  implicit none

#include "constants.h"
  
  character(len=1), intent(in)       :: char
  
  return
end subroutine Logfile_break
