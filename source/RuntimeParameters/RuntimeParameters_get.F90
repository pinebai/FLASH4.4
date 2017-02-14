!!****f* source/RuntimeParameters/RuntimeParameters_get
!!
!! NAME
!!  RuntimeParameters_get
!!
!! SYNOPSIS
!!
!!  RuntimeParameters_get(  char*(in) :: name,
!!              real/int/str/log(out) :: value) 
!!
!! DESCRIPTION
!!
!!  This function retrievess a runtime parameter from a linked list.
!!  RuntimeParameters_get is an overloaded routine.  Underneath the
!!  hood, RuntimeParameters_get implements RuntimeParameters_getReal,
!!  RuntimeParameters_getInt, RuntimeParameters_getStr and
!!  RuntimeParameters_getLog.
!!
!!  Typically RuntimeParameters_get should only be called in the
!!  initialization routines for each unit (Driver_init, IO_init,
!!  Hydro_init, Simulation_init etc).  In FLASH3 we get the runtime
!!  parameters one time and store the values in the corresponding data
!!  modules (Driver_data, IO_data, Hydro_data, Simulation_data, etc.)
!!
!! ARGUMENTS
!!
!!  name:       name of parameter
!!  value:      parameter value
!!
!!
!! EXAMPLE
!!
!!  use RuntimeParameters, ONLY : RuntimeParameters_set
!!  
!!     integer :: lrefine_max
!!
!!     call RuntimeParameters_get('lrefine_max', lrefine_max) 
!!     if (lrefine_max .GT. 2) then
!!          .......
!!
!! NOTES
!!   
!!  Because RuntimeParameters_get is an overloaded function, a user calling
!!  the routine must USE the interface RuntimeParameters_interface.
!!
!!
!!
!!***
   
subroutine RuntimeParameters_getReal (name, value)
implicit none
  character(len=*), intent(in)          :: name
  real, intent(out)                     :: value
end subroutine RuntimeParameters_getReal

subroutine RuntimeParameters_getInt (name, value)
implicit none
  character(len=*), intent(in)           :: name
  integer, intent(out)                   :: value
end subroutine RuntimeParameters_getInt

subroutine RuntimeParameters_getStr (name, value)
implicit none
  character(len=*),intent(in)             :: name
  character(len=*),intent(out)            :: value
end subroutine RuntimeParameters_getStr

subroutine RuntimeParameters_getLog (name, value)
implicit none
  character(len=*),intent(in)              :: name
  logical,intent(out)                      :: value
end subroutine RuntimeParameters_getLog

