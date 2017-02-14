!!****f* source/monitors/Profiler/Profiler_start
!!
!! NAME
!!  Profiler_start
!!
!! SYNOPSIS
!!
!!  Profiler_start(character(IN) :: name, OR
!!                 integer(IN)   :: id)
!!                   
!!  
!! DESCRIPTION 
!!  
!!  Start profiling the section 'name' or 'id'
!!  
!! ARGUMENTS 
!!
!!  name - the name of the section to profile
!!  id - the integer id of the section to profile
!!
!!***

subroutine Profiler_startName(name)
implicit none
  character (len=*), intent(in) :: name

end subroutine Profiler_startName

subroutine Profiler_startId(id)
implicit none
  integer, intent(in) :: id

end subroutine Profiler_startId
