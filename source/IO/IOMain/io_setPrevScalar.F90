!!****if* source/IO/IOMain/io_setPrevScalar
!!
!! NAME
!!  io_setPrevScalar
!!
!! SYNOPSIS
!!
!!  io_setPrevScalar(char*(in) :: name,
!!                       real/int/str/log(in) :: value) 
!!
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
!!  io_setPrevScalar is typically called by each Unit's outputScalar
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
!!  Because io_setPrevScalar is an overloaded function a user calling
!!  the routine must include the header file IO.h
!!  
!!  (Under the hood io_setPrevScalar under the hood it calls 
!!  io_setPrevScalarReal, io_setPrevScalarInt, io_setPrevScalarStr, 
!!  or io_setPrevScalarLog and keeps a separate list for each type.)
!!
!!
!! EXAMPLE
!! 
!!  To checkpoint simulation time Driver_outputScalars calls
!!
!!  #include "IO.h"
!!
!!  call io_setPrevScalar("time", simTime)
!!
!!***

subroutine io_setPrevScalarReal (name, value)

  use IO_data, only : io_scalar
  
implicit none
  character(len=*), intent(in)          :: name
  real, intent(in)                      :: value
  logical                               :: current_val = .FALSE.

  call NameValueLL_setReal(io_scalar, name, value, current_val)

  return

end subroutine io_setPrevScalarReal


  

  
subroutine io_setPrevScalarInt (name, value)
  
  use IO_data, only : io_scalar

implicit none
  character(len=*), intent(in)           :: name
  integer, intent(in)                    :: value
  logical                               :: current_val = .FALSE.

    
  call NameValueLL_setInt(io_scalar, name, value, current_val)
  
  return
  
end subroutine io_setPrevScalarInt



subroutine io_setPrevScalarStr (name, value)

  use IO_data, only : io_scalar
  
implicit none
  character(len=*),intent(in)             :: name, value
  logical                               :: current_val = .FALSE.
  
  call NameValueLL_setStr(io_scalar, name, value, current_val)

  return
  
end subroutine io_setPrevScalarStr






subroutine io_setPrevScalarLog (name, value)
  
  use IO_data, only : io_scalar

implicit none
  character(len=*),intent(in)              :: name
  logical,intent(in)                       :: value
  logical                               :: current_val = .FALSE.
  
  call NameValueLL_setLog(io_scalar, name, value, current_val)
  
  return
  
end subroutine io_setPrevScalarLog


