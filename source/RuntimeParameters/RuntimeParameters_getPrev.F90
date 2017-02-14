!!****f* source/RuntimeParameters/RuntimeParameters_getPrev
!!
!! NAME
!!  RuntimeParameters_getPrev
!!
!! SYNOPSIS
!!
!!  RuntimeParameters_getPrev(char*(in) :: name,
!!                     real/int/str/log(out) :: value) 
!!                      
!!
!! DESCRIPTION
!!
!! This function gets an initial parameter from a linked list
!! implemented under the hood.  'Previous' refers to the parameter
!! values that are used to reset the state of the simulation when
!! starting from restart.  For example, the checkpoint file might
!! store lrefine_max = 5, but in the flash.par a user restarting the
!! code may want to set lrefine_max to 6 to increase the resolution of
!! the grid.  The code must be restarted with the grid in the same
!! state as it was when it was checkpointed.  However, after that,
!! various values can be changed.  Adding the routines getPrev*
!! allows the user to do some checks to see if the checkpoint
!! parameter and the flash.par parameter values are the same.
!!
!! The getPrev* routines should only be used for these type of
!! comparisons.  If a user simply wants to get a parameter use
!! "RuntimeParameters_get"
!!
!! ARGUMENTS
!!
!! name:       name of parameter
!! value:      parameter value
!!
!!
!! NOTES
!!   
!!  Because RuntimeParameters_getPrev is an overloaded function, a user calling
!!  the routine must USE the interface RuntimeParameters_interface. 
!!
!!
!!***
   
subroutine RuntimeParameters_getPrevReal (name, value)
implicit none
  character(len=*), intent(in)          :: name
  real, intent(out)                     :: value
end subroutine RuntimeParameters_getPrevReal

subroutine RuntimeParameters_getPrevInt (name, value)
implicit none
  character(len=*), intent(in)           :: name
  integer, intent(out)                    :: value
end subroutine RuntimeParameters_getPrevInt

subroutine RuntimeParameters_getPrevStr (name, value)
implicit none
  character(len=*),intent(in)             :: name
  character(len=*),intent(out)             :: value
end subroutine RuntimeParameters_getPrevStr

subroutine RuntimeParameters_getPrevLog (name, value)
implicit none
  character(len=*),intent(in)              :: name
  logical,intent(out)                       :: value
end subroutine RuntimeParameters_getPrevLog

