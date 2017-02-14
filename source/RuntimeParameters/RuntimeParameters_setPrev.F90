!!****f* source/RuntimeParameters/RuntimeParameters_setPrev
!!
!! NAME
!!  RuntimeParameters_setPrev
!!
!! SYNOPSIS
!!
!!  RuntimeParameters_setPrev(char*(in) :: name,
!!                     real/int/str/log(in) :: value) 
!!
!! DESCRIPTION
!!
!!  This is an overloaded function. RuntimeParameters_setPrev
!!  under the hood implements RuntimeParameters_setPrevReal
!!  RuntimeParameters_setPrevInt, RuntimeParameters_setPrevStr
!!  RuntimeParameters_setPrevLog
!!
!!  This routine allows the user to keep track of a value from a restart run.
!!  The linked List holding the runtime parameters can hold 2 values, 
!!  a current value, and an "initial" or "previous" value.  The routine sets
!!  the value of the latter.  
!!
!!
!!
!! EXAMPLE
!!
!!  Take a run using Paramesh with a refinement level of 4.  The user
!!  runs the simulation for some time then writes a checkpoint. After
!!  examining the checkpoint data the user wants to increase the maximum
!!  resolution in the simulation and sets the refinement level to 6.
!!  Upon restart the grid may need to compare the values 4 and 6 and refine
!!  the grid further.  The value '4' would be read in from the checkpoint
!!  file and stored in the 'previous', while the value '6' would
!!  be stored in the 'current' slot.  We imagine there are only a few 
!!  runtime parameters in which this routine would be relevant.
!!
!!  use RuntimeParameters, ONLY : RuntimeParameters_set
!!  
!!     integer :: lrefine_max
!!
!!     lrefine_max = 4
!!     call RuntimeParameters_setPrev('lrefine_max', lrefine_max) 
!!
!!
!! ARGUMENTS
!!
!!  name:       name
!!  value:      name value
!!
!!
!!
!! NOTES
!!   
!!  Because RuntimeParameters_setPrev is an overloaded function, a user calling
!!  the routine must USE the interface RuntimeParameters_interface. 
!!
!!***

subroutine RuntimeParameters_setPrevReal (name, value)
implicit none
  character(len=*), intent(in)          :: name
  real, intent(in)                      :: value
  logical                               :: current_val = .FALSE.
end subroutine RuntimeParameters_setPrevReal
  
subroutine RuntimeParameters_setPrevInt (name, value)
implicit none
  character(len=*), intent(in)           :: name
  integer, intent(in)                    :: value
  logical                               :: current_val = .FALSE.
end subroutine RuntimeParameters_setPrevInt

subroutine RuntimeParameters_setPrevStr (name, value)
implicit none
  character(len=*),intent(in)             :: name, value
  logical                               :: current_val = .FALSE.
end subroutine RuntimeParameters_setPrevStr

subroutine RuntimeParameters_setPrevLog (name, value)
implicit none
  character(len=*),intent(in)              :: name
  logical,intent(in)                       :: value
  logical                               :: current_val = .FALSE.
end subroutine RuntimeParameters_setPrevLog


