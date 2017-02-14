!!****f* source/monitors/Profiler/Profiler_stop
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
implicit none
  character (len=*), intent(in) :: name

end subroutine Profiler_stopName

subroutine Profiler_stopId(id)
implicit none
  integer, intent(in) :: id

end subroutine Profiler_stopId
