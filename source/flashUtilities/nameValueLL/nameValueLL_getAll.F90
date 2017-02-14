!!****if* source/flashUtilities/nameValueLL/nameValueLL_getAll
!!
!! NAME
!!  nameValueLL_getAll
!!
!! SYNOPSIS
!!
!!  nameValueLL_getAll (context_type(IN)                    :: context,
!!                      integer(INOUT)                      :: num,
!!                      character(OUT)                      :: names(num),
!!                      real/integer/character/logical(OUT) :: values(num),
!!                      logical(OUT)                        :: changed(num))
!!
!! DESCRIPTION
!!
!! Returns the names and values for a given context (parameter or scalar) 
!! of a given type (int, real, str, log) and if they have changed since 
!! the start of the simulation. 
!!
!! NOTE: in case of restarting the change is checked against the values
!!       in the checkpoint file and not the DEFAULT values in config file
!!
!! from a linked list implemented under
!! the hood.  This subroutine is overloaded and implements
!! nameValueLL_getAllReal, nameValueLL_getAllInt
!! nameValueLL_getAllStr, nameValueLL_getAllLog
!!
!! ARGUMENTS
!! context:    type of list (ie for parameters or scalars)
!! num:        number of parameters or scalars to get
!! names:       names of parameters or scalars
!! values:      value
!! changed:    logicals indicating if parameter changed or not 
!!             from the initial run
!!
!!***


subroutine nameValueLL_getAllReal(context, num, names, values, changed)
  
  use nameValueLL_data !, ONLY: context_type, nameValueLL_find, &
 !     &  name_invalid, name_real, name_int, name_str, name_log, &
 !     &  real_list_type, int_list_type, str_list_type, log_list_type

  implicit none
#include "constants.h"

  type (context_type),intent(in)        :: context
  integer, intent(INOUT) :: num
  character(len=MAX_STRING_LENGTH), intent(OUT) :: names(num)
  real, intent(OUT) :: values(num)
  logical, intent(OUT) :: changed(num)
  type(real_list_type), pointer :: node
  integer :: i,max_return,num_returned
  
  max_return = num
  i = 1
  num_returned = 0
  node => context%real_list
  do while (associated(node).and.(num_returned < max_return))
     names(i) = node%name
     values(i) = node%value
     if ( node%value .eq. node%initValue ) then
        changed(i) = .false.
     else
        changed(i) = .true.
     endif
     i = i + 1
     num_returned = num_returned + 1
     node => node%next
  enddo
  num = num_returned
end subroutine nameValueLL_getAllReal


subroutine nameValueLL_getAllInt(context, num, names, values, changed)
  
  use nameValueLL_data !, ONLY: context_type, nameValueLL_find, &
 !     &  name_invalid, name_real, name_int, name_str, name_log, &
 !     &  real_list_type, int_list_type, str_list_type, log_list_type

  implicit none
#include "constants.h"

  type (context_type),intent(in)                  :: context   
  integer, intent(INOUT) :: num
  character(len=MAX_STRING_LENGTH), intent(INOUT) :: names(num)
  integer, intent(INOUT) :: values(num)
  logical, intent(OUT) :: changed(num)
  type(int_list_type), pointer :: node
  integer :: i, max_return, num_returned
  
  max_return = num
  i = 1
  num_returned = 0
  node => context%int_list
  do while (associated(node).and.(num_returned < max_return))
     names(i) = node%name
     values(i) = node%value
     if ( node%value .eq. node%initValue ) then
        changed(i) = .false.
     else
        changed(i) = .true.
     endif
     i = i + 1
     num_returned = num_returned + 1
     node => node%next
  enddo

  num = num_returned
  

end subroutine nameValueLL_getAllInt

subroutine nameValueLL_getAllStr(context, num, names, values, changed)

  use nameValueLL_data !, ONLY: context_type, nameValueLL_find, &
 !     &  name_invalid, name_real, name_int, name_str, name_log, &
  !    &  real_list_type, int_list_type, str_list_type, log_list_type

  implicit none
#include "constants.h"

  type (context_type),intent(in)                :: context   
  integer, intent(INOUT) :: num
  character(len=MAX_STRING_LENGTH), intent(OUT) :: names(num)
  character(len=MAX_STRING_LENGTH), intent(OUT) :: values(num)
  logical, intent(OUT) :: changed(num)
  type(str_list_type), pointer :: node
  integer :: i,max_return,num_returned
  
  max_return = num
  i = 1
  num_returned = 0
  node => context%str_list
  do while (associated(node).and.(num_returned < max_return))
     names(i) = node%name
     values(i) = node%value
     if ( node%value .eq. node%initValue ) then
        changed(i) = .false.
     else
        changed(i) = .true.
     endif
     i = i + 1
     num_returned = num_returned + 1
     node => node%next
  enddo
  num = num_returned
end subroutine nameValueLL_getAllStr

subroutine nameValueLL_getAllLog (context, num, names, values, changed)

  use nameValueLL_data !, ONLY: context_type, nameValueLL_find, &
 !     &  name_invalid, name_real, name_int, name_str, name_log, &
 !     &  real_list_type, int_list_type, str_list_type, log_list_type

implicit none
#include "constants.h"

  type (context_type),intent(in)                :: context   
  integer, intent(INOUT) :: num
  character(len=MAX_STRING_LENGTH), intent(OUT) :: names(num)
  logical, intent(OUT) :: values(num)
  logical, intent(OUT) :: changed(num)
  type(log_list_type), pointer :: node
  integer :: i,max_return,num_returned
  
  max_return = num
  i = 1
  num_returned = 0
  node => context%log_list
  do while (associated(node).and.(num_returned < max_return))
     names(i) = node%name
     values(i) = node%value
!! some compilers (cube) doesn't like checking for equality of logical values 
!! so we work around it
     if ( node%value .and. node%initValue ) then
        changed(i) = .false.
     else if ( node%value .or. node%initValue ) then
        changed(i) = .true.
     else
        changed(i) = .false.
     endif
     i = i + 1
     num_returned = num_returned + 1
     node => node%next
  enddo
  num = num_returned
end subroutine nameValueLL_getAllLog

