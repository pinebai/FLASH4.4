!!****if* source/monitors/Timers/TimersMain/MPINative/tmr_create
!!
!! NAME
!!  tmr_create
!!
!! SYNOPSIS
!!
!!  tmr_create(character(IN) :: name, 
!!                 integer(OUT), optional :: i)
!!
!! DESCRIPTION
!!  Create a new timer, and optionally return the
!!  integer key to that timer
!!
!! ARGUMENTS
!!  name --    a string containing the name of the section to time
!!  i --       an optional integer parameter that will return the
!!              integer key to the timer.
!!
!! PARAMETERS
!!
!!***

subroutine tmr_create (name,i)

  use Timers_data, ONLY: tmr_timerInvalid
  implicit none

! Arguments
  character(len=*), intent(IN)  :: name
  integer,optional, intent(OUT) :: i
! Local variables
  integer          :: id
  
#ifdef NOOP
  return
#endif
  
  call tmr_findTimerIndex (name,.FALSE., id)
  
  if (id == tmr_timerInvalid) then
     call tmr_findTimerIndex (name,.TRUE., id)
  endif
  
  if (present(i)) i = id
  
  return
  
end subroutine tmr_create
