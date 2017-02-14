!!****if* source/monitors/Timers/TimersMain/MPINative/Timers_start
!!
!! NAME
!!  Timers_start - start a timer given a string name
!!
!! SYNOPSIS
!!
!!  Timers_start(character(IN) :: name(:))
!!
!! DESCRIPTION
!!  Start the timer specified by a supplied name
!!
!! ARGUMENTS
!!  name --   a string containing the name of the timer to start
!!
!! PARAMETERS
!!
!! NOTES
!!
!!  Multithreaded timers
!!
!!  * Do not call Timers_start and Timers_stop within an OpenMP single
!!    section because timing information is only recorded for the
!!    master thread.  It is possible that the single thread is the
!!    master thread on some MPI ranks but a non-master thread on other
!!    MPI ranks.  When this happens the timer stacks will be different
!!    across MPI ranks and we will get the 'whaddayamakeofthat?'
!!    message.
!!
!!  * Do not replace omp master with an omp single in the Timers unit
!!    because it may cause a deadlock.
!!
!!  * Do not replace omp master with a nowait omp single in the Timers
!!    unit because the Timers unit is of unknown thread-safety.  This
!!    would be a problem if different threads are in different nowait
!!    single sections.
!!
!!***

subroutine Timers_startString(name)
use Timers_interface, ONLY : Timers_startIndex

!  use Timers_data

  implicit none

  character(len=*), intent(in) :: name
  
  integer            :: i
  
#ifdef NOOP
  return
#endif

  !$omp master
  call tmr_findTimerIndex (name,.TRUE., i)
  !$omp end master
  call Timers_startIndex(i)

  return

end subroutine Timers_startString

!!****if* source/monitors/Timers/TimersMain/Timers_start
!!
!! NAME
!!   Timers_start - start a timer given a integer key
!!
!! SYNOPSIS
!!
!!   Timers_startIndex(integer(IN) :: i)
!!
!! DESCRIPTION
!!   Start timing the timer specified by a supplied integer key
!!
!! ARGUMENTS
!!   i --     an integer key specifiying the timer to start
!!
!! PARAMETERS
!!
!!***

subroutine Timers_startIndex (i)
 
  use Timers_data, ONLY: tmr_numSegments, tmr_acctSegs, tmr_callStack
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

  integer, intent(in) :: i

  integer :: j, pushResult
  integer :: currentStackIndex
  real :: temp_time

  !$omp master
  if ((i < 1) .or. (i > tmr_numSegments)) then
     write (*,*) ' Timers_startIndex: error: invalid timer ', i
  else 
     call tmr_stackListIndex(tmr_acctSegs(i)%stacks, tmr_callStack, j)
     if (j==0) then
        ! The next comment is an sgi IPA compiler directive: don't inline the 
        ! stackList_add function.  The compiler inlined incorrectly, and 
        ! the code broke.  Can't figure out why.
        
        !*$* NOINLINE here (tmr_stackListSdd)
        call tmr_stackListAdd(tmr_acctSegs(i)%stacks, tmr_callStack, j)
     end if
     call tmr_stackPush(tmr_callStack, i, pushResult)
     if (pushResult < 0) then
        call Driver_abortFlash('[Timers_start] Ran out of space on timer call stack. ' //&
             'Probably means calling start without a corresponding stop.')
     end if
     if (j < 0) then
        write (*,*) 'perfmon: ran out of space for timer, "', trim(tmr_acctSegs(i)%name), '", cannot time this timer with perfmon'
     else
        currentStackIndex = j
        tmr_acctSegs(i)%isTimed(currentStackIndex) = .true.
        call tmr_etime(temp_time)
        tmr_acctSegs(i)%dtime(currentStackIndex) = temp_time
        tmr_acctSegs(i)%timesCalled(currentStackIndex) = tmr_acctSegs(i)%timesCalled(currentStackIndex) + 1
!!$          call Profiler_start(tmr_acctSegs(i)%name, i)
     end if
  end if
  !$omp end master
  return
end subroutine Timers_startIndex
