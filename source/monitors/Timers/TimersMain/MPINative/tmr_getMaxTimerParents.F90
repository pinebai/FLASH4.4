!!****if* source/monitors/Timers/TimersMain/MPINative/tmr_getMaxTimerParents
!!
!! NAME
!!   tmr_getMaxTimerParents
!!
!! SYNOPSIS
!!   tmr_getMaxTimerParents(integer(OUT) :: value)
!!
!! DESCRIPTION
!!  Accessor function for how many places a timer can appear 
!!  in various call stacks.  For instance, timers A, B, and C 
!!  appear in a run, as follows:
!!
!!   /A
!!   /B/A
!!   /C/B/A
!!   /B/C
!!
!!  then A has three parents, B has two, and C has two.  This number is 
!!  in general statically limited.
!!
!! ARGUMENTS
!!   value -- the number of places a timer can appear
!!
!! PARAMETERS
!!
!!***

subroutine tmr_getMaxTimerParents(inmaxParents)

  use Timers_data, ONLY: tmr_maxTimerParents

  implicit none

  integer, intent(OUT) :: inmaxParents

  inmaxParents = tmr_maxTimerParents
  return

end subroutine tmr_getMaxTimerParents
