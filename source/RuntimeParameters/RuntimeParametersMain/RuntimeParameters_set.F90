!!****if* source/RuntimeParameters/RuntimeParametersMain/RuntimeParameters_set
!!
!! NAME
!!  RuntimeParameters_set
!!
!! SYNOPSIS
!!
!!  RuntimeParameters_set(char*(in) :: name,
!!               real/int/str/log(in) :: value) 
!!
!! DESCRIPTION
!!
!!  This is an overloaded function. RuntimeParameters_set under the
!!  hood implements RuntimeParameters_setReal
!!  RuntimeParameters_setInt, RuntimeParameters_setStr
!!  RuntimeParameters_setLog.  
!!
!!  Parameters are added to the runtime database.
!!
!!  This routine is normally called from RuntimeParameters_read and
!!  other internal IO routines.  We have left it in the API interface
!!  for user flexibility.
!!  
!! ARGUMENTS
!!
!!  name:       name
!!  value:      name value
!!
!!
!! EXAMPLE
!!
!!  use RuntimeParameters, ONLY : RuntimeParameters_set
!!  
!!     integer :: lrefine_max
!!
!!     lrefine_max = 6
!!     call RuntimeParameters_set('lrefine_max', lrefine_max) 
!!
!! NOTES
!!   
!!  Because RuntimeParameters_set is an overloaded function, a user calling
!!  the routine must include the header file RuntimeParameters.h
!!
!!
!!
!!!!
!!***

subroutine RuntimeParameters_setReal (name, value)

  use RuntimeParameters_data, only : parameter
  
implicit none
  character(len=*), intent(in)          :: name
  real, intent(in)                      :: value
  logical                               :: current_val = .TRUE.

  call nameValueLL_setReal(parameter, name, value, current_val)

  return

end subroutine RuntimeParameters_setReal


  

  
subroutine RuntimeParameters_setInt (name, value)
  
  use RuntimeParameters_data, only : parameter

implicit none
  character(len=*), intent(in)           :: name
  integer, intent(in)                    :: value
  logical                               :: current_val = .TRUE.
  
  call nameValueLL_setInt(parameter, name, value, current_val)
  
  return
  
end subroutine RuntimeParameters_setInt



subroutine RuntimeParameters_setStr (name, value)

  use RuntimeParameters_data, only : parameter
  
implicit none
  character(len=*),intent(in)             :: name, value
  logical                               :: current_val = .TRUE.
  
  call nameValueLL_setStr(parameter, name, value, current_val)

  return
  
end subroutine RuntimeParameters_setStr






subroutine RuntimeParameters_setLog (name, value)
  
  use RuntimeParameters_data, only : parameter

implicit none
  character(len=*),intent(in)              :: name
  logical,intent(in)                       :: value
  logical                               :: current_val = .TRUE.

  
  call nameValueLL_setLog(parameter, name, value, current_val)
  
  return
  
end subroutine RuntimeParameters_setLog


