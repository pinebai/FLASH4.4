!!****f* source/monitors/Logfile/Logfile_close
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

  implicit none

  logical, optional, intent(IN) :: logUnitLocal

end subroutine Logfile_close

