!!****if* source/monitors/Timers/TimersMain/MPINative/tmr_findTimerIndex
!!
!! NAME
!!   tmr_findTimerIndex - find a timer with a given name
!!
!! SYNOPSIS
!!   tmr_findTimerIndex(character(:), intent(IN) :: name, 
!!                      logical, intent(IN) :: createIfNone,
!!                      integer, intent(OUT) :: index)
!!
!! DESCRIPTION
!!   Find the integer key for a given name name, optionally 
!!   create it if it doesn't exist
!!
!! ARGUMENTS
!!   name --         a string containing the name of the timer to find
!!   createIfNone --  a logical determining if the routine is to create
!!                    a timer if it can't find it.
!!   index --        index of the timer
!!
!! PARAMETERS
!!
!!***

subroutine tmr_findTimerIndex (name, createIfNone, index)
  use Timers_data, ONLY: tmr_acctSegs, tmr_numSegments, tmr_maxSegments, &
       tmr_timerInvalid

  implicit none  
  character(len=*), intent(in) :: name
  logical, intent(in)          :: createIfNone
  integer, intent(out) :: index
  integer            :: i
  
  i = 1
  do while (i <= tmr_numSegments)
     if (tmr_acctSegs(i)%name == name) exit
     i = i + 1
  enddo
  
  if ((i > tmr_numSegments) .and. (tmr_numSegments < tmr_maxSegments-1)) then
     if (createIfNone) then
!!$          call profile_define(name,i)
        tmr_numSegments = tmr_numSegments + 1
        tmr_acctSegs(i)%name     = name
        tmr_acctSegs(i)%time(1)     = 0.
        tmr_acctSegs(i)%dtime(1)    = 0.
        tmr_acctSegs(i)%pctTime(1) = 0.
        tmr_acctSegs(i)%isTimed(1) = .false.
     else
        i = tmr_timerInvalid
     endif
  endif
  index = i
  
  return
end subroutine tmr_findTimerIndex
