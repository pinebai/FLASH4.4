!!****if* source/diagnostics/ProtonImaging/localAPI/pi_openMonitorFile
!!
!! NAME
!!
!!  pi_openMonitorFile
!!
!! SYNOPSIS
!!
!!  call pi_openMonitorFile ()
!!
!! DESCRIPTION
!!
!!  Opens the monitor file for recording proton imaging progress. The file name
!!  has already been set during initialization.
!!
!! ARGUMENTS
!!
!!  none
!!
!! NOTES
!!          
!!  Only the master processor opens the file.
!!
!!***

subroutine pi_openMonitorFile ()

  implicit none

  return
end subroutine pi_openMonitorFile
