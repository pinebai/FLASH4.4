!!****if* source/monitors/Timers/TimersMain/Tau/Timers_reset
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
!!   Reset the accumulated data in timers.  This does not delete
!!   timers, so integer mappings to named timers will remain intact,
!!   but all information previously gathered about timers will be
!!   lost.  In the middle of a run, for instance, this could be called
!!   once per timestep along with Timers_getSummary to get timer
!!   summary information for each timestep.
!!  
!! ARGUMENTS 
!!
!!***

subroutine Timers_reset()
  
  use Driver_interface, ONLY : Driver_abortFlash
  implicit none

  call Driver_abortFlash &
       ("[Timers_reset]: Not appropriate in this timer implementation")

  return
end subroutine Timers_reset

