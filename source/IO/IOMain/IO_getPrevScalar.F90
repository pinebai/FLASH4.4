!!****if* source/IO/IOMain/IO_getPrevScalar
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
!! error:      error code
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
!!  An error code of 0 (NORMAL) indicates normal execution.  A code of -1 
!!  (NOTFOUND) indicates that the scalar value was not found in the checkpoint
!!  file.  This case must be handled by the calling subroutine if the optional
!!  error dummy argument was provided.
!!
!!***


#include "constants.h"
   
subroutine IO_getPrevScalarReal (name, value, error)

  use IO_data, ONLY : io_scalar
  use Driver_interface, ONLY : Driver_abortFlash

implicit none
  character(len=*), intent(in)          :: name
  real, intent(out)                      :: value
  logical                        :: current_val = .FALSE.
  integer, intent(out),optional :: error
  integer :: retError

  if (present(error)) error = 0
  
  call NameValueLL_getReal(io_scalar, name, value, current_val, retError)
  
  if (present(error)) then 
     error = retError
  else
     if(retError /= NORMAL) then
        call Driver_abortFlash("[IO_getPrevScalar] ERROR: Could not find real scalar value!")
     end if
  end if

  return
  
end subroutine IO_getPrevScalarReal


   
subroutine IO_getPrevScalarInt (name, value, error)

  use IO_data, ONLY : io_scalar
  use Driver_interface, ONLY : Driver_abortFlash

implicit none
  character(len=*), intent(in)          :: name
  integer, intent(out)                :: value
  logical                        :: current_val = .FALSE.
  integer, intent(out),optional :: error
  integer :: retError  

  if (present(error)) error = 0
  
  call NameValueLL_getInt(io_scalar, name, value, current_val, retError)
 if (present(error)) then 
     error = retError
  else
     if(retError /= NORMAL) then
        call Driver_abortFlash("[IO_getPrevScalar] ERROR: Could not find integer scalar value!")
     end if
  end if

  return
  
end subroutine IO_getPrevScalarInt


   
subroutine IO_getPrevScalarStr (name, value, error)

  use IO_data, ONLY : io_scalar
  use Driver_interface, ONLY : Driver_abortFlash

implicit none
  character(len=*), intent(in)          :: name
  character(len=*), intent(out)         :: value
  logical                        :: current_val = .FALSE.
  integer, intent(out),optional :: error
  integer :: retError

  if (present(error)) error = 0

  call NameValueLL_getStr(io_scalar, name, value, current_val, retError)
 if (present(error)) then 
     error = retError
  else
     if(retError /= NORMAL) then
        call Driver_abortFlash("[IO_getPrevScalar] ERROR: Could not find string scalar value!")
     end if
  end if

  return
  
end subroutine IO_getPrevScalarStr


   
subroutine IO_getPrevScalarLog (name, value, error)

  use IO_data, ONLY : io_scalar
  use Driver_interface, ONLY : Driver_abortFlash

implicit none
  character(len=*), intent(in)          :: name
  logical, intent(out)                      :: value
  logical                        :: current_val = .FALSE.
  integer, intent(out), optional :: error
  integer :: retError

  if (present(error)) error = 0
  
  call NameValueLL_getLog(io_scalar, name, value, current_val, retError)
 if (present(error)) then 
     error = retError
  else
     if(retError /= NORMAL) then
        call Driver_abortFlash("[IO_getPrevScalar] ERROR: Could not find logical scalar value!")
     end if
  end if

  return
  
end subroutine IO_getPrevScalarLog



