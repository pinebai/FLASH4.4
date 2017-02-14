!!****ih* source/monitors/Timers/TimersMain/Tau/tmr_interface
!!
!!***
Module tmr_interface
#include "Flash.h"
#include "constants.h"

  interface tmr_findTimerIndex
     subroutine tmr_findTimerIndex(name, createIfNone, index) 
       character(len=*), intent(in) :: name
       logical, intent(in)          :: createIfNone
       integer, intent(out) :: index
     end subroutine tmr_findTimerIndex
  end interface

  interface tmr_etime
     subroutine tmr_etime(time)
       real, intent(OUT) :: time
     end subroutine tmr_etime
  end interface

end Module tmr_interface
