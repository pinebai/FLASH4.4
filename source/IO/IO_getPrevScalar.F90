!!****f* source/IO/IO_getPrevScalar
!!
!!
!! NAME
!!  IO_getPrevScalar
!!
!! SYNOPSIS
!!
!!  IO_getPrevScalar(char*(in) :: name,
!!                   real/int/str/log(out) :: value, 
!!                   integer(out),optional :: error)
!!                      
!!
!! DESCRIPTION
!!
!!  This function gets a 'previous' scalar from a 
!!  linked list implemented under the hood.  
!!  The scalar list can store 2 values, a 'current' value and
!!  a 'previous' value.  Previous values are used when restarting 
!!  a run and might also be used to compare
!!  a scalar value stored in a checkpoint to a new value set
!!  apon restart. (see example below)
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
!!  comment, the name of the logfile, or the setup line.
!!
!!
!! ARGUMENTS
!!
!! name:       name of scalar
!! value:      scalar value
!! error:      error code. 0 (NORMAL) for normal termination
!!
!! EXAMPLE
!!
!!  When restarting a run a user might want to change the timestep, dt
!!  or adjust the simulation time. The checkpoint read routine reads the
!!  scalar values from the checkpoint file and stores them in the 'previous'
!!  slot. The user can then compare the 'previous' value of the scalar to
!!  the new value.
!!
!!     USE IO_interface : ONLY IO_getPrevScalar
!!  
!!     real :: dt
!!
!!     call IO_getPrevScalar('dt', dt) 
!!
!! NOTES
!!   
!!  Because IO_getPrevScalar is an overloaded function a user calling
!!  the routine must use IO_interface.F90.
!!  
!!  (Under the hood, a call to the generic name IO_getPrevScalar is handled
!!  by IO_getPrevScalarReal, IO_getPrevScalarInt, IO_getPrevScalarStr, 
!!  or IO_getPrevScalarLog. A separate list is kept for each type.)
!!
!!  An error value of NOTFOUND (-1) indicated that the given previous value
!!  was not found.  The calling subroutine should handle this error if the
!!  optional error dummy argument was provided.
!!  
!!***

   
subroutine IO_getPrevScalarReal (name, value, error)


implicit none
  character(len=*), intent(in)          :: name
  real, intent(out)                      :: value
  integer, intent(out),optional :: error


  if(present(error)) error = 0

  value = 0.0
  return
  
end subroutine IO_getPrevScalarReal


   
subroutine IO_getPrevScalarInt (name, value, error)


implicit none
  character(len=*), intent(in)          :: name
  integer, intent(out)                :: value
  integer, intent(out), optional :: error

  if(present(error)) error = 0

  value = 0
  return
  
end subroutine IO_getPrevScalarInt


   
subroutine IO_getPrevScalarStr (name, value, error)


implicit none
  character(len=*), intent(in)          :: name
  character(len=*), intent(out)         :: value
  integer, intent(out),optional :: error

  if(present(error)) error = 0

  value = ' '
  return
  
end subroutine IO_getPrevScalarStr


   
subroutine IO_getPrevScalarLog (name, value, error)


implicit none
  character(len=*), intent(in)          :: name
  logical, intent(out)                      :: value
  integer, intent(out), optional :: error

  if(present(error)) error = 0

  value = .FALSE.
  return
  
end subroutine IO_getPrevScalarLog



