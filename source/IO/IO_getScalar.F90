!!****f* source/IO/IO_getScalar
!!
!! NAME
!!  IO_getScalar
!!
!! SYNOPSIS
!!
!!  IO_getScalar(char*(in) :: name,
!!               real/int/str/log(in) :: value) 
!!                      
!!
!! DESCRIPTION
!!
!!  This function gets a scalar from a 
!!  linked list implemented under the hood.
!!   
!!  In FLASH3 what we mean by scalars are single value variables
!!  associated with an entire flash run.
!!  These scalars are in contrast to Grid scope variables which
!!  need to be stored at each zone of each block in the simulation.
!!  Density, pressure and temperature are examples of Grid scope 
!!  variables while simTime, dt, globalNumBlocks are single quantities
!!  associated with the entire run. 
!!
!!  Scalars of this type can be integers, reals, strings or logical
!!  values.  An example of a string scalar might be the FLASH3 run
!!  comment, name of the logfile or setup line
!!
!! ARGUMENTS
!!
!! name:       name of scalar
!! value:      scalar value
!!
!!
!!
!! EXAMPLE
!!
!!     USE IO_interface, ONLY : IO_getScalar
!!  
!!     real :: dt
!!
!!     call IO_getScalar('dt', dt) 
!!
!! NOTES
!!   
!!  Because IO_getScalar is an overloaded function a user calling
!!  the routine must USE the interface IO_interface.
!!  
!!  (Under the hood IO_getScalar under the hood it calls 
!!  IO_getScalarReal, IO_getPreviousScalarInt, IO_getScalarStr, 
!!  or IO_getScalarLog and keeps a separate list for each type.)
!!
!!
!!***


   
subroutine IO_getScalarReal (name, value)


implicit none
  character(len=*), intent(in)          :: name
  real, intent(inout)                      :: value

  return
  
end subroutine IO_getScalarReal


   
subroutine IO_getScalarInt (name, value)


implicit none
  character(len=*), intent(in)          :: name
  integer, intent(inout)                :: value

  return
  
end subroutine IO_getScalarInt


   
subroutine IO_getScalarStr (name, value)


implicit none
  character(len=*), intent(in)          :: name, value
  return
  
end subroutine IO_getScalarStr


   
subroutine IO_getScalarLog (name, value)


implicit none
  character(len=*), intent(in)          :: name
  logical, intent(inout)                      :: value

  return
  
end subroutine IO_getScalarLog



