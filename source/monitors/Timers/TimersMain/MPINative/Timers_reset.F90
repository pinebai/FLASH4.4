!!****if* source/monitors/Timers/TimersMain/MPINative/Timers_reset
!!
!! NAME
!!  Timers_reset
!!
!! SYNOPSIS
!!
!!  Timers_reset()
!!  
!! DESCRIPTION 
!!  
!!  Reset the accumulated data in timers.  This does not delete
!!  timers, so integer mappings to named timers will remain intact,
!!  but all information previously gathered about timers will be
!!  lost.  In the middle of a run, for instance, this could be called
!!  once per timestep along with Timers_getSummary to get timer
!!  summary information for each timestep.
!!  
!! ARGUMENTS 
!!
!!***

subroutine Timers_reset()

  use Timers_data, ONLY: tmr_acctSegs, tmr_callStack, tmr_maxSegments, &
       tmr_maxTimerParents, tmr_numSegments, tmr_initDate, tmr_initTime
  implicit none

  integer   :: i, j

  !This initialization is here so that TAU generates valid code.
  !Without it, TAU declares a variable after !$omp master.
  i = 0

  !$omp master
  call tmr_stackZero(tmr_callStack)
  do i = 1, tmr_maxSegments

     call tmr_stackListZero(tmr_acctSegs(i)%stacks)
     do j = 1, tmr_maxTimerParents
        tmr_acctSegs(i)%time(j)     = 0.
        tmr_acctSegs(i)%dtime(j)    = 0.
        tmr_acctSegs(i)%timesCalled(j) = 0
        tmr_acctSegs(i)%isTimed(j) = .false.
     enddo
  enddo

  
  call current_date_time(tmr_initDate)
  call tmr_etime(tmr_initTime)
  !$omp end master
!!$  ! initialize the profiling library
!!$    call profile_initialize()

end subroutine Timers_reset

