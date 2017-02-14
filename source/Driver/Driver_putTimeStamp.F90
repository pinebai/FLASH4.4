!!****f* source/Driver/Driver_putTimeStamp
!!
!! NAME
!!   Driver_putTimeStamp
!! 
!! SYNOPSIS
!!
!!   Driver_putTimeStamp(character(40), intent(IN) :: dateStr)
!!
!! DESCRIPTION
!!
!!   Sets the string which contains the time stamp of the run.
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


subroutine Driver_putTimeStamp(dateStr)

  implicit none
  character(len=40), intent(IN)     :: dateStr

  return

end subroutine Driver_putTimeStamp
