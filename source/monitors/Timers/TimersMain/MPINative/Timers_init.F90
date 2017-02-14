!!****if* source/monitors/Timers/TimersMain/MPINative/Timers_init
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
!!
!!***

#include "constants.h"
subroutine Timers_init( initialWCTime)

  use Timers_data, ONLY: tmr_initDate, tmr_initTime, tmr_writeStatSummary, &
       tmr_eachProcWritesSummary, tmr_globalMe, tmr_globalNumProcs, tmr_globalComm
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_getMype, Driver_getNumProcs, Driver_getComm

  implicit none

  real, intent(out) :: initialWCTime

  !This initialization is here so that TAU generates valid code.
  !Without it, TAU declares a variable after !$omp master.
  initialWCTime = 0.0

  ! Everybody should know this
  !$omp master
  call Driver_getComm(GLOBAL_COMM, tmr_globalComm)
  call Driver_getMype(GLOBAL_COMM, tmr_globalMe)
  call Driver_getNumProcs(GLOBAL_COMM, tmr_globalNumProcs)

  call RuntimeParameters_get("writeStatSummary", tmr_writeStatSummary)
  call RuntimeParameters_get("eachProcWritesSummary", tmr_eachProcWritesSummary)

  call tmr_etime(initialWCTime)
  tmr_initTime = initialWCTime
  call current_date_time(tmr_initDate)

  call tmr_init()
  !$omp end master

end subroutine Timers_init
