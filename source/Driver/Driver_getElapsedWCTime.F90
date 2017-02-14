!!****f* source/Driver/Driver_getElapsedWCTime
!!
!! NAME
!!  Driver_getElapsedWCTime
!!
!! SYNOPSIS
!!  
!!
!!  Driver_getElapsedWCTime(real(out) :: elapsedWCTime)
!!  
!! DESCRIPTION 
!!
!!  This is an accessor funtion that returns the elapsed wall clock time
!!  since the beginning of the run.
!!
!! ARGUMENTS
!!
!!  elapsedWCTime - returned value of elapsed wall clock time in seconds
!!
!!
!!***


!! NOTES
!!
!! The Driver unit uses a few unit scope variables that are
!! accessible to all routines within the unit, but not to the
!! routines outside the unit. These variables begin with "dr_"
!! like, dr_globalMe or dr_dt, dr_beginStep, and are stored in FORTRAN
!! module Driver_data (in file Driver_data.F90). The other variables
!! are local to the specific routine and do not have the prefix "dr_"

subroutine Driver_getElapsedWCTime(elapsedWCTime)

implicit none
  real, intent(out) :: elapsedWCTime

  !dummy values for stubs
  elapsedWCTime = 0.0
  
end subroutine Driver_getElapsedWCTime

