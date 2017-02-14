!!****f* source/monitors/Timers/Timers_start
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
  implicit none

  character(len=*), intent(in) :: name

  return
end subroutine Timers_startString

!!****if* source/monitors/Timers/Timers_startIndex
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
  implicit none
  integer, intent(in) :: i
  return
end subroutine Timers_startIndex
