!!****f* source/IO/IO_setScalar
!!
!! NAME
!!  IO_setScalar
!!
!! SYNOPSIS
!!
!!  IO_setScalar(char*(in) :: name,
!!               real/int/str/log(in) :: value) 
!!
!! DESCRIPTION
!!
!!  Accessor routine that sets a scalar value to a scalar list which
!!  will then be checkpointed or written to a plotfile.
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
!!  IO_setScalar is typically called by each Unit's outputScalar
!!  (ie Driver_outputScalars, Grid_outputScalars) 
!!  routines right before checkpointing or writing a plotfile.
!!  A user wishing to write a new scalar to a checkpoint file would
!!  need to call this routine
!!
!! ARGUMENTS
!!
!!  name:       name
!!  value:      name value
!!
!! 
!! NOTES
!!   
!!  Because IO_setScalar is an overloaded function a user calling
!!  the routine must USE the interface IO_interface.
!!  
!!  (Under the hood, IO_setScalar calls IO_setScalarReal,
!!  IO_setScalarInt, IO_setScalarStr, or IO_setScalarLog and keeps a
!!  separate list for each type.)
!!
!!
!! EXAMPLE
!! 
!!  To checkpoint the simulation time, Driver_sendOutputData does
!!  call IO_setScalar("time", dr_simTime)
!!
!!***

subroutine IO_setScalarReal (name, value)

implicit none
  character(len=*), intent(in)          :: name
  real, intent(in)                      :: value

  return

end subroutine IO_setScalarReal


  

  
subroutine IO_setScalarInt (name, value)
  
implicit none
  character(len=*), intent(in)           :: name
  integer, intent(in)                    :: value
  
  return
  
end subroutine IO_setScalarInt



subroutine IO_setScalarStr (name, value)
  
implicit none
  character(len=*),intent(in)             :: name, value

  return
  
end subroutine IO_setScalarStr






subroutine IO_setScalarLog (name, value)

implicit none
  character(len=*),intent(in)              :: name
  logical,intent(in)                       :: value
  
  return
  
end subroutine IO_setScalarLog


