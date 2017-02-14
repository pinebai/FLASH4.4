!!****if* source/monitors/Timers/TimersMain/MPINative/tmr_getMaxCallStackDepth
!!
!! NAME
!!   tmr_getMaxCallStackDepth
!!
!! SYNOPSIS
!!   tmr_getMaxCallStackDepth(integer(OUT) :: value)
!!
!! DESCRIPTION
!!  Accessor function for how deep the timers can nest at runtime, which is 
!!  in general, a fixed depth.
!!
!! ARGUMENTS
!!   value -- the max depth for nesting timers.
!!
!! PARAMETERS
!!
!!***

subroutine tmr_getMaxCallStackDepth(value)

  use Timers_data, ONLY: tmr_maxCallStackDepth

  implicit none

  integer, intent(OUT) :: value

  value = tmr_maxCallStackDepth
  return

end subroutine tmr_getMaxCallStackDepth
