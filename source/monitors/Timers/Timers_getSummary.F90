!!****f* source/monitors/Timers/Timers_getSummary
!!
!! NAME
!!   
!! Timers_getSummary
!!
!! SYNOPSIS
!!
!!  call Timers_getSummary(integer(in) :: nIntervals)
!!
!! DESCRIPTION
!!  Write out the Timers summary to the logfile.
!!  Timers collects the performance data and then calls
!!  Logfile_writeSummary to do the actual formatting.
!!
!!  It should be safe to call this at any point in the 
!!  simulation-- an updated summary will appear in the 
!!  logfile.  
!!
!! ARGUMENTS
!!
!!  nIntervals - number of subintervals timed, which is determined 
!!               by the caller of this code, but this will 
!!               typically be the number of timesteps taken
!!               since the last time Timers_init was called. 
!!
!!***


subroutine Timers_getSummary( nIntervals)

  implicit none

  integer, intent(in) :: nIntervals

  return

end subroutine Timers_getSummary



