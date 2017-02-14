!!****f* source/monitors/Profiler/Profiler_getSummary
!! NAME
!!   
!! Profiler_getSummary
!!
!! SYNOPSIS
!!  Profiler_getSummary(integer(in) :: nIntervals)
!!
!! DESCRIPTION
!!   Write out the Profiler summary.  This can be done in whatever
!!   format the particular profiler chooses to do it.  
!!
!! ARGUMENTS
!!
!!  nIntervals - nIntervals    :: number of subintervals timed
!!                                 typically timsteps
!! 
!!
!!
!!***


subroutine Profiler_getSummary(nIntervals)

  implicit none
  integer, intent(in) :: nIntervals

end subroutine Profiler_getSummary

