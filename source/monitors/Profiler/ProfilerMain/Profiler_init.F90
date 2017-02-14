!!****if* source/monitors/Profiler/ProfilerMain/Profiler_init
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

subroutine Profiler_init()
  use Profiler_data, ONLY : prf_evolutionOnly
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  implicit none
  call RuntimeParameters_get("profileEvolutionOnly", prf_evolutionOnly)
end subroutine Profiler_init
