!!****if* source/monitors/Timers/TimersMain/MPINative/Timers_stop
!!
!! NAME
!!  Timers_stopIndex - stop a timer given an integer key
!!
!! SYNOPSIS
!!
!!  Timers_stopIndex(integer(IN) :: i)
!!
!! DESCRIPTION
!!  Stop timing a timer specified by a supplied integer key.  This
!!  implementation also calls Profiler_stop with the corresponding
!!  timer name.
!!
!! ARGUMENTS
!!  i --   an integer key specifiying the timer to stop
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

recursive subroutine Timers_stopIndex (i)

  use Timers_data, ONLY:  tmr_callStack, tmr_acctSegs
  use Timers_interface, ONLY : Timers_startIndex

  implicit none  

  integer, intent(in) :: i
  integer            :: result
  integer            :: j
  real :: temp_time

#ifdef NOOP
  return
#endif
  !$omp master
  call tmr_stackTop(tmr_callStack, result)
  if (i == result) then
     call tmr_stackPop(tmr_callStack, j)
     if (j >= 0) then
!!$       call profile_end(tmr_acctSegs(i)%name, i)      
        call tmr_stackListIndex(tmr_acctSegs(i)%stacks, tmr_callStack, j)
        call tmr_etime(temp_time)
        tmr_acctSegs(i)%time(j) = temp_time - tmr_acctSegs(i)%dtime(j) + tmr_acctSegs(i)%time(j)
        tmr_acctSegs(i)%isTimed(j) = .false.
     end if
  else
     call tmr_stackTop(tmr_callStack, j)
     if (j >= 0) then
        call Timers_stopIndex(j)
        call Timers_stopIndex(i)
        call Timers_startIndex(j)
     end if
  end if
  !$omp end master
  return

end subroutine Timers_stopIndex

!!****if* source/Timers/TimersMain/Timers_stopString
!!
!! NAME
!!   Timers_stopString
!!
!! SYNOPSIS
!!
!!   Timers_stopString(character(IN) :: name(:))
!!
!! DESCRIPTION
!!
!!   Stop timing a timer specified by a supplied name. This
!!   implementation also calls Profiler_stop.
!!
!! ARGUMENTS
!!   name --   a name specifiying the timer to stop
!!
!! PARAMETERS
!!
!!***
subroutine Timers_stopString(name)
  use Timers_interface, ONLY : Timers_stopIndex

  implicit none  

  character(len=*), intent(IN)   :: name
  integer            :: i
  
  !$omp master
  call tmr_findTimerIndex (name,.FALSE., i)
  !$omp end master
  call Timers_stopIndex(i)
  
  return
end subroutine Timers_stopString
