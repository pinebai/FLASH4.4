!!****if* source/monitors/Timers/TimersMain/MPINative/tmr_lookupIndex
!!
!! NAME
!!   tmr_lookupIndex 
!!
!! SYNOPSIS
!!
!!   call tmr_lookupIndex( character(len=*), intent(IN)  :: name, 
!!                              int, intent(OUT)         :: index)
!!
!! DESCRIPTION
!!   find the integer key (for faster use) for a given name
!!
!! ARGUMENTS
!!
!!   name --  a string containing the name of a timer
!!   index -- the returned timer index
!!
!! PARAMETERS
!!
!!***

subroutine tmr_lookupIndex(name, index)

  implicit none  

  character(len=*), intent(IN) :: name
  integer, intent(out) :: index
  
  call tmr_findTimerIndex(name,.FALSE., index)
  
  return
  
end subroutine tmr_lookupIndex
