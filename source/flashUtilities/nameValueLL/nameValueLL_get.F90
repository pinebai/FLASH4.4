!!****if* source/flashUtilities/nameValueLL/nameValueLL_get
!!
!! NAME
!!  nameValueLL_get
!!
!! SYNOPSIS
!!
!!  nameValueLL_get (context_type(IN)                      :: context,
!!                   character(IN)                         :: name,
!!                   real/integer/character/logical(INOUT) :: value,
!!                   logical(IN)                           :: current_val)
!!
!! DESCRIPTION
!!
!! Gets a name from a linked list implemented under
!! the hood.  This subroutine is overloaded and implements
!! nameValueLL_getReal, nameValueLL_getInt
!! nameValueLL_getStr, nameValueLL_getLog
!!
!! ARGUMENTS
!! context:    type of list (ie for parameters or scalars)
!! name:       name or parameter or scalar
!! value:      name value
!! current_val: a logical true/false value that indicates whether you
!!          want the current value of the parameter or scalar or if
!!          you want the original value (ie the value that was stored
!!          in the checkpoint file)
!!
!! NOTES
!!    For most purposes current_val will be set to true.  This argument is
!!    hidden under the hood from users.  Users will call either
!!    RuntimeParameters_get or RuntimeParameters_getPrev.  The rare times
!!    that users call RuntimeParameters_getPrev, the logical value current_value
!!    will be set to false.
!!
!!***

   
subroutine nameValueLL_getReal (context, name, value, current_val, error)

  use nameValueLL_data !, ONLY: context_type, nameValueLL_find, &
 !     &  name_invalid, name_real, name_int, name_str, name_log, &
 !     &  real_list_type, int_list_type, str_list_type, log_list_type
  use Driver_interface, ONLY : Driver_abortFlash
  use Logfile_interface, ONLY : Logfile_stampMessage

  implicit none

#include "constants.h"
#include "Flash_mpi.h"

  type (context_type),intent(in)        :: context
  character(len=*), intent(in)          :: name
  real, intent(inout)                   :: value
  logical, intent(in)                   :: current_val
  type (real_list_type), pointer    :: node
  character(len=MAX_STRING_LENGTH)  :: buff1, buff2
  integer, intent(out) :: error
  integer :: myPE, ierr

  call MPI_comm_rank(MPI_COMM_WORLD, myPE, ierr)

  error = 0

  call nameValueLL_find (context, name, node)
  
  if (associated(node)) then
     if (current_val) then
        value = node%value
     else if(.NOT. current_val) then
        value = node%initValue
     else
        
        call Driver_abortFlash("nameValueLL_getReal: logical value current_val has improper value")
     end if
  else
     buff1 = "WARNING: requested real parameter '" // trim(name) // "' not found"
     call Logfile_stampMessage(buff1)
     if (myPE == MASTER_PE) print *,buff1
     error = NOTFOUND
     !call Driver_abortFlash(buff1)
  endif
  
  return
  
end subroutine nameValueLL_getReal



subroutine nameValueLL_getInt (context, name, value, current_val, error)

  use nameValueLL_data !, ONLY: context_type, nameValueLL_find, &
 !     &  name_invalid, name_real, name_int, name_str, name_log, &
 !     &  real_list_type, int_list_type, str_list_type, log_list_type
  use Driver_interface, ONLY : Driver_abortFlash
  use Logfile_interface, ONLY : Logfile_stampMessage

  implicit none
#include "constants.h"
#include "Flash_mpi.h"

  type (context_type),intent(in)        :: context
  character(len=*), intent(in)           :: name
  integer, intent(inout)                 :: value
  logical, intent(in)                    :: current_val
  type (int_list_type), pointer          :: node
  character(len=MAX_STRING_LENGTH)  :: buff1, buff2
  integer, intent(out)  :: error

  integer :: myPE, ierr
  
  call MPI_comm_rank(MPI_COMM_WORLD, myPE, ierr)
  
  error = NORMAL

  call nameValueLL_find (context, name, node)
  if (associated(node)) then
     if (current_val) then
        value = node%value
     else if(.NOT. current_val) then
        value = node%initValue
     else
        call Driver_abortFlash("nameValueLL_getInt: logical value current_val has improper value")
     end if
  else
     buff1 = "WARNING: requested integer parameter '" // trim(name) // "' not found"
     call Logfile_stampMessage(buff1)
     if(myPE == MASTER_PE) print *, buff1
     error = NOTFOUND
     !call Driver_abortFlash(buff1)
  endif
  
  return
  
end subroutine nameValueLL_getInt


subroutine nameValueLL_getStr (context, name, value, current_val, error)

  use nameValueLL_data !, ONLY: context_type, nameValueLL_find, &
 !     &  name_invalid, name_real, name_int, name_str, name_log, &
 !     &  real_list_type, int_list_type, str_list_type, log_list_type
  use Driver_interface, ONLY : Driver_abortFlash
  use Logfile_interface, ONLY : Logfile_stampMessage

  implicit none
#include "constants.h"
#include "Flash_mpi.h"

  type (context_type),intent(in)        :: context
  character(len=*),intent(in)             :: name
  character(len=*),intent(inout)          :: value
  logical, intent(in)                     :: current_val
  type (str_list_type), pointer           :: node
  character(len=MAX_STRING_LENGTH)  :: buff1, buff2
  integer, intent(out) :: error
  integer :: myPE, ierr

  call MPI_comm_rank(MPI_COMM_WORLD, myPE, ierr)

  error = NORMAL
  
  call nameValueLL_find (context, name, node)
  if (associated(node)) then
     if (current_val) then
        value = node%value
     else if(.NOT. current_val) then
        value = node%initValue
     else
        call Driver_abortFlash("nameValueLL_getStr: logical value current_val has improper value")
     end if
  else
     buff1 = "WARNING: requested string parameter '" // trim(name) // "' not found"
     call Logfile_stampMessage(buff1)
     if(myPE == MASTER_PE) print *, buff1
     error = NOTFOUND
     !call Driver_abortFlash(buff1)
  endif
  
  return
  
end subroutine nameValueLL_getStr






subroutine nameValueLL_getLog (context, name, value, current_val, error)

  use nameValueLL_data !, ONLY: context_type, nameValueLL_find, &
 !     &  name_invalid, name_real, name_int, name_str, name_log, &
 !     &  real_list_type, int_list_type, str_list_type, log_list_type
  use Driver_interface, ONLY : Driver_abortFlash
  use Logfile_interface, ONLY : Logfile_stampMessage

  implicit none
#include "constants.h"
#include "Flash_mpi.h"

  type (context_type),intent(in)        :: context
  character(len=*),intent(in)            :: name
  logical,intent(inout)                  :: value
  logical, intent(in)                    :: current_val
  type (log_list_type), pointer          :: node
  character(len=MAX_STRING_LENGTH)  :: buff1, buff2
  integer, intent(out)  :: error
  integer :: myPE, ierr

  call MPI_comm_rank(MPI_COMM_WORLD, myPE, ierr)

  error = 0
  

  call nameValueLL_find (context, name, node)
  if (associated(node)) then
     if (current_val) then
        value = node%value
     else if(.NOT. current_val) then
        value = node%initValue
     else
        call Driver_abortFlash("nameValueLL_getLog: logical value current_val has improper value")
     end if
  else
     buff1 = "WARNING: requested logical parameter '" // trim(name) // "' not found"
     call Logfile_stampMessage(buff1)
     if(myPE == MASTER_PE) print *, buff1
     error = NOTFOUND
     !call Driver_abortFlash(buff1)
  endif
  
  return
end subroutine nameValueLL_getLog



