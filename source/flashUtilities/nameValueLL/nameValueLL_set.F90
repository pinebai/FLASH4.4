!!****if* source/flashUtilities/nameValueLL/nameValueLL_set
!!
!! NAME
!!  nameValueLL_set
!!
!! SYNOPSIS
!!
!!  nameValueLL_set(context_type(INOUT)  :: context,
!!                  character(IN)        :: name,
!!                  real/int/str/log(IN) :: value,
!!                  logical(IN)          :: current_val)
!!
!! DESCRIPTION
!!
!! Sets real parameters into a linked list.  This is an overloaded
!! function nameValueLL_set under the hood implements
!! nameValueLL_setReal nameValueLL_setInt, nameValueLL_setStr
!! nameValueLL_setLog
!!
!! ARGUMENTS
!!
!! context:      structure holding all the data
!! name:         name to be set
!! value:        value of name to be set
!! current_val:  distinguishes between restart (false) and initial (true) values
!!
!!***

subroutine nameValueLL_setReal (context, name, value, current_val)

  use nameValueLL_data !, ONLY: context_type, &
 !     &  nameValueLL_find, nameValueLL_check, nameValueLL_addReal, &
 !     &  name_invalid, name_real, name_int, name_str, name_log, &
 !     &  real_list_type, int_list_type, str_list_type, log_list_type, TYPE_VAR
  use Driver_interface, ONLY : Driver_abortFlash
  use Logfile_interface, ONLY : Logfile_stamp
  
#include "constants.h"
implicit none


  type (context_type),intent(inout)                   :: context
  character(len=*), intent(in)          :: name
  real, intent(in)                      :: value
  logical, intent(in)                   :: current_val
  type (real_list_type),pointer     :: node
  logical                           :: valid
  character(len=MAX_STRING_LENGTH)       :: buf
  character(len=300) :: logStr
  character(len=80) :: initValueStr, valueStr
  

  call nameValueLL_find (context, name, node)
  
  if (associated(node)) then
     if (.NOT.  node%isConstant) then
        call nameValueLL_check(node,value,valid)
        if (.not. valid) then
           call nameValueLL_logRulesReal(name,node%numValues,node%minValues,node%maxValues)
           write (buf,'(F8.3)') value
           call Driver_abortFlash("nameValue_set: Trying to set '"// trim(name) //"' to invalid value "// trim(buf))
        endif
        if (current_val) then
           node%value = value
        else
           node%initValue = value
        endif
     else if (current_val .AND. (node%value.EQ.value)) then
        return    ! Same value as existing value marked as constant - RETURN
     else if (.NOT. current_val) then
        if(node%initValue.NE.value) then
           write(initValueStr, '(es20.13)') node%initValue
           write(valueStr, '(es20.13)') value
           write(logStr, "(a,a,a,a)") &
                'current=', trim(adjustl(initValueStr)), &
                ', new(checkpoint)=', trim(adjustl(valueStr))
           call Logfile_stamp(logStr, &
                '[nameValueLL_set] Different previous value for '//trim(name))
           node%initValue = value
        end if
     else
        write(*,*) "set : Can not change parameter with constant attribute:", name
        call Driver_abortFlash('ERROR: unable to change constant parameter')
     endif
  else
     ! name is not found - add it to list 
     call nameValueLL_addReal(context, name, value, TYPE_VAR)
  endif
  
  return    
end subroutine nameValueLL_setReal


  

  
subroutine nameValueLL_setInt (context, name, value, current_val)
  
  use nameValueLL_data !, ONLY: context_type, &
 !     &  nameValueLL_find, nameValueLL_check, nameValueLL_addInt, &
 !     &  name_invalid, name_real, name_int, name_str, name_log, &
 !      &  real_list_type, int_list_type, str_list_type, log_list_type, TYPE_VAR
  use Driver_interface, ONLY : Driver_abortFlash
  use Logfile_interface, ONLY : Logfile_stamp

#include "constants.h"

implicit none

  type (context_type), intent(inout)       :: context             
  character(len=*), intent(in)           :: name
  integer, intent(in)                    :: value
  logical, intent(in)                   :: current_val
  type (int_list_type),pointer       :: node
  logical                              :: valid
  character(len=MAX_STRING_LENGTH)     :: buf
  character(len=300) :: logStr
  character(len=80) :: initValueStr, valueStr

  call nameValueLL_find (context, name, node)
  
  if (associated(node)) then
     if (.NOT. node%isConstant) then
        call nameValueLL_check(node,value,valid)
        if (.not. valid) then
           call nameValueLL_logRulesInt(name,node%numValues,node%minValues,node%maxValues)
           write (buf,'(I12)') value
           call Driver_abortFlash("Trying to set '"// trim(name) // "' to invalid value " // trim(buf))
        endif
        if (current_val) then
           node%value = value
        else
           node%initValue = value
        endif
     else if (current_val .AND. (node%value.EQ.value)) then
        return    ! Same value as existing value marked as constant - RETURN
     else if (.NOT. current_val) then
        if(node%initValue.NE.value) then
           write(initValueStr, '(i12)') node%initValue
           write(valueStr, '(i12)') value
           write(logStr, "(a,a,a,a)") &
                'current=', trim(adjustl(initValueStr)), &
                ', new(checkpoint)=', trim(adjustl(valueStr))
           call Logfile_stamp(logStr, &
                '[nameValueLL_set] Different previous value for '//trim(name))
           node%initValue = value
        end if
     else
        write(*,*) "set : Can not change name with constant attribute:", name
        call Driver_abortFlash('ERROR: unable to change constant name')
     end if
  else
     ! could not find name so we are adding it to list
     call nameValueLL_addInt(context, name, value, TYPE_VAR)
  endif
  
  return
  
end subroutine nameValueLL_setInt




subroutine nameValueLL_setStr (context, name, value, current_val)
  
  use nameValueLL_data !, ONLY: context_type, &
 !     &  nameValueLL_find, nameValueLL_check, nameValueLL_addStr, &
 !     &  name_invalid, name_real, name_int, name_str, name_log, &
 !     &  real_list_type, int_list_type, str_list_type, log_list_type, TYPE_VAR
  use Driver_interface, ONLY : Driver_abortFlash
  use Logfile_interface, ONLY : Logfile_stamp

#include "constants.h"

implicit none

  type (context_type), intent(inout)       :: context            
  character(len=*),intent(in)             :: name, value
  logical, intent(in)                   :: current_val
  type (str_list_type),pointer        :: node
  logical                              :: valid
  character(len=300) :: logStr

  call nameValueLL_find (context, name, node)
  
  if (associated(node)) then
     if (.NOT. node%isConstant) then
        call nameValueLL_check(node,value,valid)
        if (.not. valid) then
           call nameValueLL_logRulesStr(name,node%numValues,node%validValues)
           call Driver_abortFlash("nameValue_set: Trying to set '"// trim(name) // "' to invalid value '"// trim(value) //"'")
        endif
        if (current_val) then
           node%value = value
        else
           node%initValue = value
        endif
     else if (current_val .AND. (node%value.EQ.value)) then
        return    ! Same value as existing value marked as constant - RETURN
     else if (.NOT. current_val) then
        if(node%initValue.NE.value) then
           write(logStr, "(a,a,a,a)") &
                'current=', trim(adjustl(node%initValue)), &
                ', new(checkpoint)=', trim(adjustl(value))
           call Logfile_stamp(logStr, &
                '[nameValueLL_set] Different previous value for '//trim(name))
           node%initValue = value
        end if
     else
        write(*,*) "set : Can not change name with constant attribute:", name
        call Driver_abortFlash('ERROR: unable to change constant name')
     end if
  else
     ! could not find name so we are adding it to list 
     call nameValueLL_addStr(context, name, value, TYPE_VAR)
  endif
    
  return
  
end subroutine nameValueLL_setStr






subroutine nameValueLL_setLog (context, name, value, current_val)
  
  use nameValueLL_data !, ONLY: context_type, &
 !     &  nameValueLL_find, nameValueLL_check, nameValueLL_addLog, &
 !     &  name_invalid, name_real, name_int, name_str, name_log, &
 !     &  real_list_type, int_list_type, str_list_type, log_list_type, TYPE_VAR
  use Driver_interface, ONLY : Driver_abortFlash
  use Logfile_interface, ONLY : Logfile_stamp

#include "constants.h"

implicit none
  type (context_type), intent(inout)       :: context            
  character(len=*),intent(in)              :: name
  logical,intent(in)                       :: value
  logical, intent(in)                   :: current_val
  type (log_list_type),pointer         :: node
  character(len=300) :: logStr
  character(len=80) :: initValueStr, valueStr

  call nameValueLL_find (context, name, node)
    
  if (associated(node)) then
     if (.NOT. node%isConstant) then
        if (current_val) then
           node%value = value
        else
           node%initValue = value
        endif
     else if (current_val .AND. (node%value.EQV.value)) then
        return    ! Same value as existing value marked as constant - RETURN
     else if (.NOT. current_val) then
        if(node%initValue.NEQV.value) then
           write(initValueStr, '(l12)') node%initValue
           write(valueStr, '(l12)') value
           write(logStr, "(a,a,a,a)") &
                'current=', trim(adjustl(initValueStr)), &
                ', new(checkpoint)=', trim(adjustl(valueStr))
           call Logfile_stamp(logStr, &
                '[nameValueLL_set] Different previous value for '//trim(name))
           node%initValue = value
        end if
     else
        write(*,*) "set: Can not change name with constant attribute:", name
        call Driver_abortFlash('ERROR: unable to change constant name')
     end if
  else
     ! could not find name so we are adding it to list 
     call nameValueLL_addLog(context, name, value, TYPE_VAR)
  endif
  
  return
  
end subroutine nameValueLL_setLog


