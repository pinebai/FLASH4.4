!!****if* source/monitors/Profiler/ProfilerMain/Profiler_data
!!
!! NAME
!!  Profiler_data
!!
!! SYNOPSIS
!!
!!  use Profiler_data
!!
!! DESCRIPTION
!!
!!  Holds the data needed by the Profiler unit
!!
!!***

module Profiler_data
  implicit none
  character (len=*), parameter :: prf_evolutionName = 'FLASH_evolution'
  logical, save :: prf_evolutionOnly
end module Profiler_data
