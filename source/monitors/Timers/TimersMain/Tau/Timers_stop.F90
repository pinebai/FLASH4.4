!!****if* source/monitors/Timers/TimersMain/Tau/Timers_stop
!!
!! NAME
!!   Timers_stopIndex - stop a timer given an integer key
!!
!! SYNOPSIS
!!
!!   Timers_stopIndex(integer(IN) :: i)
!!
!! DESCRIPTION
!!   Stop timing a timer specified by a supplied integer key
!!
!! ARGUMENTS
!!   i --   an integer key specifiying the timer to stop
!!
!! PARAMETERS
!!
!!***

#include "constants.h"

subroutine Timers_stopIndex (i)

  use Driver_interface, ONLY : Driver_abortFlash
  use Timers_data, ONLY : tmr_freeSlot, tmr_tauList

  implicit none  
  integer, intent(in) :: i

  !i >= tmr_freeSlot because we do not create the timer on a stop call.
  if (i < 1 .or. i >= tmr_freeSlot) then
     call Driver_abortFlash("[Timer_stopIndex]: Timer index not valid.")
  end if


  !Check the state of timerStarted.  A .false. should never happen.
  if (tmr_tauList(i) % timerStarted .eqv. .false.) then
     call Driver_abortFlash("[Timers_stopIndex]: Timer never started.")
  end if


  call TAU_PROFILE_STOP(tmr_tauList(i) % tauSavedData)
  tmr_tauList(i) % timerStarted = .false.
  
end subroutine Timers_stopIndex

!!****if* source/monitors/Timers/Tau/Timers_stopString
!!
!! NAME
!!
!!   Timers_stopString
!!
!! SYNOPSIS
!!   Timers_stopString(character(IN) :: name(len=*))
!!
!! DESCRIPTION
!!   Stop timing a timer specified by a supplied name
!!
!! ARGUMENTS
!!   name --   a name specifiying the timer to stop
!!
!! PARAMETERS
!!
!!***
subroutine Timers_stopString(name)

  use Timers_interface, ONLY : Timers_stopIndex
  use tmr_interface, ONLY : tmr_findTimerIndex
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none  

  character(len=*), intent(IN)   :: name
  integer :: index

  call tmr_findTimerIndex(name, .false., index)

  if (index == NONEXISTENT) then
     call Driver_abortFlash("[Timer_stopString]: Unable to find timer " // name)
  end if

  call Timers_stopIndex(index)

end subroutine Timers_stopString
