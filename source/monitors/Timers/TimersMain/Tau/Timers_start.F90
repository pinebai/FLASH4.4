!!****if* source/monitors/Timers/TimersMain/Tau/Timers_start
!!
!! NAME
!!   Timers_start - start a timer given a string name
!!
!! SYNOPSIS
!!
!!   Timers_start(character(IN) :: name)
!!
!! DESCRIPTION
!!   Start the timer specified by a supplied name
!!
!! ARGUMENTS
!!   name --   a string containing the name of the timer to start
!!
!! PARAMETERS
!!
!!***

subroutine Timers_startString(name)

  use Timers_interface, ONLY : Timers_startIndex
  use tmr_interface, ONLY : tmr_findTimerIndex
  implicit none

  character(len=*), intent(in) :: name
  integer :: index

  !If we do not find the input string, the .true. argument
  !adds the input string to our timer data structure.
  call tmr_findTimerIndex(name, .true., index)

  call Timers_startIndex(index)

  return
end subroutine Timers_startString

!!****if* source/monitors/Timers/Tau/Timers_startIndex
!!
!! NAME
!!   Timers_start - start a timer given a integer key
!!
!! SYNOPSIS
!!
!!   Timers_start(integer(IN) :: i)
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

  use Driver_interface, ONLY : Driver_abortFlash
  use Timers_data, ONLY : tmr_freeSlot, tmr_tauList, tmr_globalMe

  implicit none
  integer, intent(in) :: i

  if ( (i < 1) .or. (i > tmr_freeSlot) ) then
     call Driver_abortFlash("[Timer_startIndex]: Timer index not valid!")
  end if


  !Protect the user in case they issue 2 successive Timer_start() 
  !calls with no intermediate Timer_stop().
  if (tmr_tauList(i) % timerStarted .eqv. .true.) then
     print *, "[Timers_startIndex]: WARNING - Processor:", tmr_globalMe, &
          "issued 2 Timers_start() calls without a Timers_stop()"
     
     return  !RETURN because timer already running.
  end if
  

  call TAU_PROFILE_START(tmr_tauList(i) % tauSavedData)
  tmr_tauList(i) % timerStarted = .true.

  return
end subroutine Timers_startIndex
