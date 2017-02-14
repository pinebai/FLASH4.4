!!****if* source/monitors/Profiler/ProfilerMain/mpihpm/Profiler_init
!!
!! NAME
!!  Profiler_init
!!
!! SYNOPSIS
!!
!!  Profiler_init()
!!                   
!!  
!! DESCRIPTION 
!!  
!!  Initialize the profiler unit
!!  
!! ARGUMENTS 
!!
!!***

#include "Profiler.h"

subroutine Profiler_init()
  use Profiler_data, ONLY : prf_evolutionOnly
  use pr_interface, ONLY : pr_prof_control
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  implicit none
  call RuntimeParameters_get("profileEvolutionOnly", prf_evolutionOnly)
  if (prf_evolutionOnly) then
     !If we wish to profile evolution and are using gprof then
     !we must disable gprof measurements first and then re-enable measurements
     !later when we are in the evolution section of FLASH.
     call pr_prof_control(PRF_DISABLE_PROFILER)
  end if
end subroutine Profiler_init
