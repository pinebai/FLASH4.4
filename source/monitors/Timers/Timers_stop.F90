!!****f* source/monitors/Timers/Timers_stop
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

recursive subroutine Timers_stopIndex (i)

  implicit none  
  integer, intent(in) :: i
  return
end subroutine Timers_stopIndex

!!****if* source/monitors/Timers/Timers_stopString
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

  implicit none  

  character(len=*), intent(IN)   :: name
  return
end subroutine Timers_stopString
