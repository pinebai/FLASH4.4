!!****if* source/Driver/DriverMain/Driver_getTimeStamp
!!
!! NAME
!!   Driver_getTimeStamp
!! 
!! SYNOPSIS
!!
!!   Driver_getTimeStamp(character(40), intent(OUT) :: dateStr)
!!
!! DESCRIPTION
!!
!!   Get the string which contains the time stamp of the run.
!!   This variable is initialized in Driver_init, and updated
!!   in Logfile_create.  If the Logfile unit is included in a
!!   Simulation, then the time stamps in the flash.log and the
!!   flash.dat will match.  
!!     
!!
!! ARGUMENTS
!!     
!!     dateStr -- a string containing the simulation's time stamp
!!
!! NOTES
!!
!!    To output dateStr, use write(*,*) dateStr(1:len_trim(dateStr))
!!
!!***


subroutine Driver_getTimeStamp(dateStr)

  use Driver_data, ONLY:  dr_timeStamp

  implicit none
  character(len=40), intent(OUT)     :: dateStr

  dateStr = dr_timeStamp
  return

end subroutine Driver_getTimeStamp
