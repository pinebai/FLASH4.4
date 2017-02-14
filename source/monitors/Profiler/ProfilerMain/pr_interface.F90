!!****ih* source/monitors/Profiler/ProfilerMain/pr_interface
!!
!! NAME
!!  pr_interface
!!
!! SYNOPSIS
!!  use pr_interface
!!
!! DESCRIPTION
!! This is the interface module for the profiler unit.
!! 
!!***

module pr_interface
  interface
     subroutine pr_prof_control(mode) &
          bind(c,name='pr_prof_control')
       use iso_c_binding, ONLY : c_int
       integer(c_int), value, intent(in) :: mode
     end subroutine pr_prof_control
  end interface

  interface
     subroutine hpctoolkit_sampling_start() &
          bind(c,name='hpctoolkit_sampling_start')
     end subroutine hpctoolkit_sampling_start
  end interface

  interface
     subroutine hpctoolkit_sampling_stop() &
          bind(c,name='hpctoolkit_sampling_stop')
     end subroutine hpctoolkit_sampling_stop
  end interface
end module pr_interface
