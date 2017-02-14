!!****h* source/monitors/Profiler/Profiler_interface
!!
!! This is the header file for the Profiler module that defines its
!! public interfaces.  Of course, there's no implementation of Profiler yet,
!! so calling these public interfaces doesn't do much good....
!!***

Module Profiler_interface

  interface Profiler_getSummary
     subroutine Profiler_getSummary(nIntervals)
       implicit none
       integer, intent(in) :: nIntervals
     end subroutine Profiler_getSummary
  end interface

  interface Profiler_init
     subroutine Profiler_init()
       implicit none
     end subroutine Profiler_init
  end interface

  interface Profiler_start
     subroutine Profiler_startName(name)
       implicit none
       character (len=*), intent(in) :: name
     end subroutine Profiler_startName

     subroutine Profiler_startId(id)
       implicit none
       integer, intent(in) :: id
     end subroutine Profiler_startId
  end interface

  interface Profiler_stop
     subroutine Profiler_stopName(name)
       implicit none
       character (len=*), intent(in) :: name
     end subroutine Profiler_stopName

     subroutine Profiler_stopId(id)
       implicit none
       integer, intent(in) :: id
     end subroutine Profiler_stopId
  end interface

end Module Profiler_interface
