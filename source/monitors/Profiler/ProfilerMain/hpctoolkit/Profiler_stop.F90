!!****if* source/monitors/Profiler/ProfilerMain/hpctoolkit/Profiler_stop
!!
!! NAME
!!  Profiler_stop
!!
!! SYNOPSIS
!!
!!  Profiler_stop( character(IN) :: name, OR
!!                 integer(IN)   :: id)
!!                   
!!  
!! DESCRIPTION 
!!  
!!  Stop profiling the section 'name' or 'id'
!!  
!! ARGUMENTS 
!!
!!  name - the name of the section to profile
!!  id - the integer id of the section to profile
!!
!!***

subroutine Profiler_stopName(name)
  use Profiler_data, ONLY : prf_evolutionName, prf_evolutionOnly
  use pr_interface, ONLY : hpctoolkit_sampling_stop
  implicit none
  character (len=*), intent(in) :: name
  if (trim(name) == prf_evolutionName .and. prf_evolutionOnly) then
     call hpctoolkit_sampling_stop()
  end if
end subroutine Profiler_stopName

subroutine Profiler_stopId(id)
  use Driver_interface, ONLY : Driver_abortFlash
  implicit none
  integer, intent(in) :: id
  call Driver_abortFlash("Not yet implemented")
end subroutine Profiler_stopId
