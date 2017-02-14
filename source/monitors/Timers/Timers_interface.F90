!!****h* source/monitors/Timers/Timers_interface
!!
!! This is the header file for the timers module that defines its
!! public interfaces.
!!***
Module Timers_interface

  interface Timers_finalize
    subroutine Timers_finalize()
    end subroutine Timers_finalize
  end interface

  interface 
    subroutine Timers_getSummary( nIntervals)
      integer, intent(in) :: nIntervals
    end subroutine Timers_getSummary
  end interface

  interface Timers_init
    subroutine Timers_init( initialWCTime)
      real, intent(out) :: initialWCTime
    end subroutine Timers_init
  end interface

  interface Timers_reset
    subroutine Timers_reset()
    end subroutine Timers_reset
  end interface
  
  interface Timers_start
    subroutine Timers_startString(name)
      character(len=*), intent(in) :: name
    end subroutine Timers_startString

    subroutine Timers_startIndex(i)
      integer, intent(in) :: i
    end subroutine Timers_startIndex
  end interface

  interface Timers_stop
    subroutine Timers_stopString(name)
      character(len=*), intent(in) :: name   
    end subroutine Timers_stopString

    recursive subroutine Timers_stopIndex(i)
      integer, intent(in) :: i
    end subroutine Timers_stopIndex
  end interface

end Module Timers_interface
