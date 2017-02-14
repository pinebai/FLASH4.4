!!****f* source/monitors/Timers/Timers_init
!!
!! NAME
!!  Timers_init
!!
!! SYNOPSIS
!!
!!  Timers_init(real(OUT) :: initialWCTime)
!!  
!! DESCRIPTION 
!!  
!!  Initialize the timer data structures.  This will
!!  essentially delete all information previously gathered by all timers
!!  and make it safe to start timers from scratch.  In the middle of
!!  a run, for instance, this could be called once per timestep along with
!!  Timers_getSummary to get timer summary information for each timestep.
!!  
!! ARGUMENTS 
!!
!!  initialWCTime -- the initial wall clock time when this was called. 
!! 
!!***

subroutine Timers_init( initialWCTime)

  implicit none
  include "mpif.h"
  real, intent(out) :: initialWCTime

  initialWCTime = MPI_WTime()

end subroutine Timers_init
