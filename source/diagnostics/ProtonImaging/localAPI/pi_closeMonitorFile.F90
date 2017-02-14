!!****if* source/diagnostics/ProtonImaging/localAPI/pi_closeMonitorFile
!!
!! NAME
!!
!!  pi_closeMonitorFile
!!
!! SYNOPSIS
!!
!!  call pi_closeMonitorFile ()
!!
!! DESCRIPTION
!!
!!  Closes the monitor file for recording proton imaging progress.
!!
!! ARGUMENTS
!!
!!  none
!!
!! NOTES
!!          
!!  Only the master processor closes the file.
!!
!!***

subroutine pi_closeMonitorFile ()

  implicit none

  return
end subroutine pi_closeMonitorFile
