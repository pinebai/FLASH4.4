!!****ih* source/flashUtilities/nameValueLL/nameValueLL_data
!!
!!  NAME
!!    nameValueLL_data
!!
!!  SYNOPSIS
!!    
!!    nameValueLL_data is a fortran module 
!!    In addition to holding runtime data this module holds a 
!!    few overloaded functions as well.  In order for overloading
!!    to work correctly and portably in fortran the functions 
!!    need to be inside of a fortran module
!!
!!  DESCRIPTION
!!
!!    Utility module that creates a linked list to hold runtime parameters and
!!    scalars
!!
!!
!!  Arguments
!!    none
!!
!!***

module nameValueLL_data
  
#include "constants.h" 
  
 integer, parameter :: NUM_MAX_RULES = 10

 type real_list_type
     character(len=MAX_STRING_LENGTH)     :: name
     real                                 :: value
     real                                 :: initValue
     logical                              :: isConstant
     integer                              :: numValues
     real, dimension(:),pointer           :: minValues
     real, dimension(:),pointer           :: maxValues
     type (real_list_type), pointer :: next
  end type real_list_type
  
  type int_list_type
     character(len=MAX_STRING_LENGTH)   :: name
     integer                            :: value
     integer                            :: initValue
     logical                            :: isConstant
     integer                            :: numValues
     integer, dimension(:),pointer      :: minValues
     integer, dimension(:),pointer      :: maxValues
     type (int_list_type), pointer  :: next
  end type int_list_type
  
  type str_list_type
     character(len=MAX_STRING_LENGTH)   :: name
     character(len=MAX_STRING_LENGTH)   :: value
     character(len=MAX_STRING_LENGTH)   :: initValue
     logical                            :: isConstant
     integer                            :: numValues
     character(len=MAX_STRING_LENGTH), &
     & dimension(:),pointer             :: validValues
     type (str_list_type), pointer  :: next
  end type str_list_type
  
  type log_list_type
     character(len=MAX_STRING_LENGTH)   :: name
     logical                            :: value
     logical                            :: initValue
     logical                            :: isConstant
     type (log_list_type), pointer  :: next
  end type log_list_type
  
  type context_type
     type (real_list_type), pointer :: real_list => NULL()
     type (int_list_type), pointer  :: int_list  => NULL()
     type (str_list_type), pointer  :: str_list  => NULL()
     type (log_list_type), pointer  :: log_list  => NULL()
     integer                :: n_real=0, n_int=0, n_str=0, n_log=0
  end type context_type



  ! Parameter type constants, returned by inquiry routines.
  
  integer, parameter :: name_real = 1, name_int = 2, & 
       &                name_str  = 3, name_log = 4, & 
       &                name_invalid = 0
  
  integer, parameter :: TYPE_CONST = 0, TYPE_VAR = 1   
  
  integer,save       :: num_params = 0



  interface nameValueLL_initContext
     module procedure initContext
  end interface

  interface nameValueLL_getType
     module procedure getType
  end interface

  interface nameValueLL_find
     module procedure nameValueLL_findReal
     module procedure nameValueLL_findInt
     module procedure nameValueLL_findStr
     module procedure nameValueLL_findLog
  end interface
  
  interface nameValueLL_add
     module procedure nameValueLL_addReal
     module procedure nameValueLL_addInt
     module procedure nameValueLL_addStr
     module procedure nameValueLL_addLog
  end interface
      
  interface nameValueLL_check
    module procedure nameValueLL_checkReal
    module procedure nameValueLL_checkInt
    module procedure nameValueLL_checkStr
  end interface

contains
  
  !!****if* source/flashUtilities/nameValueLL/initContext
  !!
  !! NAME
  !!   initContext
  !!
  !! SYNOPSIS
  !!   call initContext(context_type(INOUT) :: context)
  !!
  !! DESCRIPTION: 
  !!    Make sure the context is in a defined initial state.
  !!
  !! ARGUMENTS
  !!    context:         name of context
  !!
  !!***
  
  subroutine initContext (context)
    implicit none
    type (context_type), intent(out)       :: context

    context%real_list => NULL()
    context%int_list  => NULL()
    context%str_list  => NULL()
    context%log_list  => NULL()
    context%n_real=0
    context%n_int=0
    context%n_str=0
    context%n_log=0
   end subroutine initContext

  !!****if* source/flashUtilities/nameValueLL/getType
  !!
  !! NAME
  !!   getType
  !!
  !! SYNOPSIS
  !!   getType(context_type(INOUT) :: context,
  !!           character(IN)       :: name,
  !!           integer(OUT)        :: name_type)
  !!
  !! DESCRIPTION: 
  !!    Given the name in a linked list and a context to search, find
  !!    the type (real, integer, etc.) of the parameter.  If the
  !!    parameter is not found, return a value which indicates that
  !!    the parameter is invalid.
  !!
  !! ARGUMENTS
  !!    context:         name of context
  !!    name:          name (in)
  !!    name_type:     type of parameter (out)
  !!
  !!***
  
  subroutine getType (context, name, name_type)
    implicit none
    type (context_type), intent(inout)       :: context
    character(len=*), intent(in)          :: name
    integer, intent(out)                  :: name_type
    type (real_list_type), pointer:: real_test
    type (int_list_type), pointer :: int_test
    type (str_list_type), pointer :: str_test
    type (log_list_type), pointer :: log_test


    call nameValueLL_find (context, name, real_test)
    call nameValueLL_find (context, name, int_test)
    call nameValueLL_find (context, name, str_test)
    call nameValueLL_find (context, name, log_test)

    name_type = name_invalid
    !               Only one of the following should test true if the filters
    !               in add_*_
    if (associated(real_test)) name_type = name_real
    if (associated(int_test))  name_type = name_int
    if (associated(str_test))  name_type = name_str
    if (associated(log_test))  name_type = name_log

    return
  end subroutine getType

  !!****if* source/flashUtilities/nameValueLL/nameValueLL_find 
  !!
  !!  NAME     
  !!   nameValueLL_find
  !!
  !!  SYNOPSIS
  !!   nameValueLL_find(context_type(IN)                        :: context,
  !!                    character(IN)                              :: name,
  !!                    real_/int_/str_/log_list_type(IN)(pointer) :: result) 
  !!
  !!  DESCRIPTION 
  !!    Find a name in a linked list using its name.
  !!    Return a pointer to the name's node if it's found,
  !!    otherwise return a null pointer. This subroutine is
  !!    overloaded it implements nameValueLL_findReal,
  !!    nameValueLL_findInt, nameValueLL_findStr,
  !!    and nameValueLL_findLog
  !!
  !!  ARGUMENTS
  !!    name:    name of string (in)
  !!    value:   name value
  !!
  !!***
  
  subroutine nameValueLL_findReal(context, name, result)
    
    type (context_type), intent(in)             :: context
    type (real_list_type), pointer              :: result
    character(len=*),intent(in)                 :: name
    character(len=len(name))                    :: name_lcase

    name_lcase = name

    call makeLowercase (name_lcase)
    result => context%real_list
    do while (associated(result))
       if (result%name == name_lcase) then
          exit
       end if
       result => result%next
    enddo
    
    return
  end subroutine nameValueLL_findReal
  
  
  subroutine nameValueLL_findInt (context, name, result)

    type (context_type), intent(in)             :: context
    type (int_list_type), pointer               :: result
    character(len=*), intent(in)                    :: name
    character(len=len(name))                        :: name_lcase
    
    name_lcase = name
    call makeLowercase (name_lcase)
    result => context%int_list
    do while (associated(result))
       if (result%name == name_lcase) exit
       result => result%next
    enddo
    
    return
  end subroutine nameValueLL_findInt
  
  
  subroutine nameValueLL_findStr(context, name, result)

    type (context_type), intent(in)             :: context    
    type (str_list_type), pointer               :: result
    character(len=*), intent(in)                    :: name
    character(len=len(name))                        :: name_lcase

    name_lcase = name
    call makeLowercase (name_lcase)
    result => context%str_list
    do while (associated(result))
       if (result%name == name_lcase) exit
       result => result%next
    enddo
    
    return
  end subroutine nameValueLL_findStr

  
  subroutine nameValueLL_findLog (context, name, result)
    
    type (context_type), intent(in)             :: context    
    type (log_list_type), pointer           :: result
    character(len=*), intent(in)                :: name
    character(len=len(name))                    :: name_lcase


    name_lcase = name
    call makeLowercase (name_lcase)
    result => context%log_list
    do while (associated(result))
       if (result%name == name_lcase) exit
       result => result%next
    enddo
    
    return
  end subroutine nameValueLL_findLog


  !!****if* source/flashUtilities/nameValueLL/nameValueLL_add
  !!
  !! NAME
  !!  nameValueLL_add
  !!
  !! SYNOPSIS
  !!   nameValueLL_add(context_type(IN)                         :: context,
  !!                   character(IN)                            :: name,
  !!                   real/integer/character/logical(pointer)  :: value) 
  !!
  !! DESCRIPTION
  !!
  !! Adds parameters to a linked list context implemented
  !!  under the hood.  This is an overloaded subroutine.
  !! It implements nameValueLL_addReal,
  !! nameValueLL_addInt, nameValueLL_addStr,
  !! nameValueLL_addLog
  !!
  !! ARGUMENTS
  !!
  !! name:       name
  !! value:      name value
  !!
  !!***
   
  subroutine nameValueLL_addReal (context, name, value, state)
    
    type (context_type), intent(inout)    :: context
    character(len=*), intent(in)          :: name
    real, intent(in)                      :: value
    integer,intent(in)                    :: state
    integer                               :: istat, name_type
    type (real_list_type), pointer    :: node, this
    
    ! Check to make sure the name doesn't already exist.

      
    call getType (context, name, name_type)
    if (name_type /= name_invalid) then
       !write (*,*) 'add :  already exists:  ', name
    else
       
       !               If it doesn't, create a node and add it to the appropriate
       !               list.
       
       allocate (node, stat=istat)
       if (istat /= 0) then
          write (*,*) 'add :  could not allocate'
       else
          node%name  = name
          node%value = value
          node%initValue = value
          node%isConstant = .false.
          node%numValues = 0
          nullify(node%minValues)
          nullify(node%maxValues)
          !! DEV: we removed the state argument -- why?
          !! reintroducing it 2009-06-05 - KW
          
          !!if (present(state)) then
             if (state == TYPE_CONST) then 
                node%isConstant = .true.
             end if
          !!end if
          call makeLowercase (node%name)
          nullify (node%next)
       endif
       if (.not. associated(context%real_list)) then
          context%real_list => node
       else
          this => context%real_list
          do while (associated(this%next))
             this => this%next
          enddo
          this%next => node
       endif
       context%n_real = context%n_real + 1
    endif
    
    return
  
  end subroutine nameValueLL_addReal
  

  subroutine nameValueLL_addInt (context, name, value, state)
    
    type (context_type), intent(inout)       :: context        
    character(len=*), intent(in)           :: name
    integer, intent(in)                    :: value
    integer, intent(in)                    :: state
    integer                                :: istat, name_type
    type (int_list_type), pointer      :: node, this
    
    !               Check to make sure the parameter doesn't already exist.
    

    
    call getType (context, name, name_type)
    if (name_type /= name_invalid) then
       !write (*,*) 'add :  already exists:  ', name
    else
       
       !            If it doesn't, create a node and add it to the appropriate
       !            list.
       
       allocate (node, stat=istat)
       if (istat /= 0) then
          write (*,*) 'add :  could not allocate'
       else
          nullify (node%next)
          node%name  = name
          node%value = value
          node%initValue = value
          node%isConstant = .false.
          node%numValues = 0
          nullify(node%minValues)
          nullify(node%maxValues)
          !!if (present(state)) then
          if (state == TYPE_CONST) then 
             node%isConstant = .true.
          end if
          !!end if
          call makeLowercase (node%name)
       endif
       if (.not. associated(context%int_list)) then
          context%int_list => node
       else
          this => context%int_list
          do while (associated(this%next))
             this => this%next
          enddo
          this%next => node
       endif
       context%n_int = context%n_int + 1
       
    endif
    return
    
  end subroutine nameValueLL_addInt
  

  subroutine nameValueLL_addStr (context, name, value, state)
    
    type (context_type), intent(inout)       :: context    
    character(len=*),intent(in)             :: name, value
    integer,intent(in)                      :: state
    integer                                 :: istat, name_type
    type (str_list_type), pointer       :: node, this
    

    ! Check to make sure the name doesn't already exist.
   
    call getType (context, name, name_type)
    if (name_type /= name_invalid) then
       !write (*,*) 'add_:  already exists:  ', name
    else
     
       !               If it doesn't, create a node and add it to the appropriate
       !               list.
       
       allocate (node, stat=istat)
       if (istat /= 0) then
          write (*,*) 'add :  could not allocate'
       else
          node%name  = name
          node%value = value
          node%initValue = value
          node%isConstant = .false.
          node%numValues = 0
          nullify(node%validValues)
          !!if (present(state)) then
             if (state == TYPE_CONST) then 
                node%isConstant = .true.
             end if
          !!end if
          call makeLowercase (node%name)
          nullify (node%next)
       endif
       if (.not. associated(context%str_list)) then
          context%str_list => node
       else
          this => context%str_list
          do while (associated(this%next))
             this => this%next
          enddo
          this%next => node
       endif
       context%n_str = context%n_str + 1
       
    endif
    
    return
    
  end subroutine nameValueLL_addStr



  subroutine nameValueLL_addLog (context, name, value, state)
    
    type (context_type), intent(inout)       :: context            
    character(len=*),intent(in)              :: name
    logical,intent(in)                       :: value
    integer,intent(in)                       :: state
    integer                                  :: istat, name_type
    type (log_list_type), pointer        :: node, this
    

    !  Check to make sure the name doesn't already exist.


    call getType (context, name, name_type)
    if (name_type /= name_invalid) then
       !write (*,*) 'add :  already exists:  ', name
    else
       
       !               If it doesn't, create a node and add it to the appropriate
       !               list.
       
       allocate (node, stat=istat)
       if (istat /= 0) then
          write (*,*) 'add :  could not allocate'
       else
          node%name  = name
          node%value = value
          node%initValue = value
          node%isConstant = .false.
          !!if (present(state)) then
             if (state == TYPE_CONST) then 
                node%isConstant = .true.
             end if
          !!end if
          call makeLowercase (node%name)
          nullify (node%next)
       endif
       if (.not. associated(context%log_list)) then
          context%log_list => node
       else
          this => context%log_list
          do while (associated(this%next))
             this => this%next
          enddo
          this%next => node
       endif
       context%n_log = context%n_log + 1
       
    endif
    
    return
    
    
  end subroutine nameValueLL_addLog


!!****if* source/flashUtilities/nameValueLL/nameValueLL_checkReal
!!
!! NAME
!!  nameValueLL_checkReal
!!
!! SYNOPSIS
!!
!!  nameValueLL_checkReal(type(real_list_type),pointer    :: node,
!!                       real                             :: value,
!!                       logical(OUT)                     :: valid)
!!
!! DESCRIPTION
!!
!! Checks if the given VALUE is valid according to rules in NODE
!!
!! ARGUMENTS
!!
!! node:       node of real_list_type
!! value:      value to be checked
!! valid:      true if valid value
!!
!!***

subroutine nameValueLL_checkReal (node, value, valid)
use Driver_interface, ONLY : Driver_abortFlash
  
implicit none

  type (real_list_type), pointer  :: node
  real, intent(in)                :: value
  logical, intent(out)            :: valid
  integer                         :: ctr
  real, parameter                 :: epsilon = TINY(1.0)
  
  valid = .false.
  if (associated(node)) then
     if (node%numValues == 0) then
        valid = .true.
     else if (.NOT.  node%isConstant) then
          do ctr = 1, node%numValues
             if ((node%minValues(ctr) .le. value +epsilon) .and. (value-epsilon .le. node%maxValues(ctr))) then
                valid = .true.
             endif
          end do
     endif
  else
     !! name is not found - add it to list 
     call Driver_abortFlash("nameValue_checkReal: invalid node given")
  endif
  
  return    
end subroutine nameValueLL_checkReal

!!****if* source/flashUtilities/nameValueLL/nameValueLL_checkInt
!!
!! NAME
!!  nameValueLL_checkInt
!!
!! SYNOPSIS
!!
!!  nameValueLL_checkInt(type(int_list_type),intent(in) :: node,
!!                       integer                        :: value,
!!                       logical(OUT)                   :: valid)
!!
!! DESCRIPTION
!!
!! Checks if the given VALUE is valid according to rules in NODE
!!
!! ARGUMENTS
!!
!! node:       node of str_list_type
!! value:      value to be checked
!! valid:      true if valid value
!!
!!***

subroutine nameValueLL_checkInt (node, value, valid)
use Driver_interface, ONLY : Driver_abortFlash

implicit none

  type (int_list_type),pointer  :: node
  integer, intent(in)             :: value
  logical, intent(out)            :: valid
  integer                         :: ctr
  
  valid = .false.
  if (associated(node)) then
     if (node%numValues == 0) then
        valid = .true.
     else if (.NOT.  node%isConstant) then
          do ctr = 1, node%numValues
             if ( (node%minValues(ctr) .le. value) .and. (value .le. node%maxValues(ctr)) ) then
                valid = .true.
             endif
          end do
     endif
  else
     !! name is not found - add it to list 
     call Driver_abortFlash("nameValue_checkAdd: invalid node given")
  endif
  
  return    
end subroutine nameValueLL_checkInt


!!****if* source/flashUtilities/nameValueLL/nameValueLL_checkStr
!!
!! NAME
!!  nameValueLL_checkStr
!!
!! SYNOPSIS
!!
!!  nameValueLL_checkStr(type(str_list_type),intent(in) :: node,
!!                       character(*),intent(in)        :: value,
!!                       logical(OUT)                   :: valid)
!!
!! DESCRIPTION
!!
!! Checks if the given VALUE is valid according to rules in NODE
!!
!! ARGUMENTS
!!
!! node:       node of str_list_type
!! value:      value to be checked
!! valid:      true if valid value
!!
!!***

subroutine nameValueLL_checkStr (node, value, valid)
use Driver_interface, ONLY : Driver_abortFlash

#include "constants.h"

implicit none

  type (str_list_type),pointer               :: node
  character(len=*), intent(in) :: value
  character(len=MAX_STRING_LENGTH)             :: lcase
  logical, intent(out)                         :: valid
  integer                                      :: ctr
  
  valid = .false.
  lcase = value
  call makeLowercase(lcase)
  if (associated(node)) then
     if (node%numValues == 0) then
        valid = .true.
     else if (.NOT.  node%isConstant) then
          do ctr = 1, node%numValues
             if (trim(node%validValues(ctr)) == trim(lcase)) then
                valid = .true.
             endif
          end do
     endif
  else
     !! name is not found - add it to list 
     call Driver_abortFlash("nameValue_checkStr: invalid node given")
  endif
  
  return    
end subroutine nameValueLL_checkStr

end module nameValueLL_data
