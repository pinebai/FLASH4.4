!!****if* source/IO/IOMain/IO_getScalar
!!
!! NAME
!!  IO_getScalar
!!
!! SYNOPSIS
!!
!!  IO_getScalar(char*(in) :: name,
!!               real/int/str/log(inout) :: value) 
!!                      
!!
!! DESCRIPTION
!!
!!  This function gets a scalar from a 
!!  linked list implemented under
!!  the hood.  
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
!!  #include "IO.h"
!!  
!!     real :: dt
!!
!!     call IO_getScalar('dt', dt) 
!!
!! NOTES
!!   
!!  Because IO_getScalar is an overloaded function a user calling
!!  the routine must include the header file IO.h
!!  
!!  (Under the hood IO_getScalar under the hood it calls 
!!  IO_getScalarReal, IO_getPreviousScalarInt, IO_getScalarStr, 
!!  or IO_getScalarLog and keeps a separate list for each type.)
!!
!!
!!***


   
subroutine IO_getScalarReal (name, value)

  use IO_data, ONLY : io_scalar
  use Driver_interface, ONLY : Driver_abortFlash

implicit none
#include "constants.h"

  character(len=*), intent(in)          :: name
  real, intent(out)                     :: value
  logical                               :: current_val = .TRUE.
  integer                               :: error



  call NameValueLL_getReal(io_scalar, name, value, current_val, error)
  if(error /= NORMAL) then
     call Driver_abortFlash("ERROR: cannot find real scalar value.")
  end if

  return
  
end subroutine IO_getScalarReal


   
subroutine IO_getScalarInt (name, value)

  use IO_data, ONLY : io_scalar
  use Driver_interface, ONLY : Driver_abortFlash

implicit none
#include "constants.h"

  character(len=*), intent(in)          :: name
  integer, intent(out)                  :: value
  logical                               :: current_val = .TRUE.
  integer                               :: error
  
  call NameValueLL_getInt(io_scalar, name, value, current_val, error)
  if(error /= NORMAL) then
     call Driver_abortFlash("ERROR: cannot find integer scalar value.")
  end if

  return
  
end subroutine IO_getScalarInt


   
subroutine IO_getScalarStr (name, value)

  use IO_data, ONLY : io_scalar
  use Driver_interface, ONLY : Driver_abortFlash

implicit none
#include "constants.h"

  character(len=*), intent(in)          :: name
  character(len=*), intent(out)         :: value
  logical                               :: current_val = .TRUE.
  integer                               :: error

  call NameValueLL_getStr(io_scalar, name, value, current_val, error)
  if(error /= NORMAL) then
     call Driver_abortFlash("ERROR: cannot find string scalar value.")
  end if

  return
  
end subroutine IO_getScalarStr


   
subroutine IO_getScalarLog (name, value)

  use IO_data, ONLY : io_scalar
  use Driver_interface, ONLY : Driver_abortFlash

implicit none
#include "constants.h"

  character(len=*), intent(in)          :: name
  logical, intent(out)                  :: value
  logical                               :: current_val = .TRUE.
  integer                               :: error

  call NameValueLL_getLog(io_scalar, name, value, current_val, error)
  if(error /= NORMAL) then
     call Driver_abortFlash("ERROR: cannot find logical scalar value.")
  end if

  return
  
end subroutine IO_getScalarLog



