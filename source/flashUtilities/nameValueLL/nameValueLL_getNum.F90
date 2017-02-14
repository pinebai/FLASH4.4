!!****if* source/flashUtilities/nameValueLL/nameValueLL_getNum
!!
!! NAME
!!  nameValueLL_getNum
!!
!! SYNOPSIS
!!
!!  nameValueLL_getNum (context_type(IN) :: context,
!!                      integer(OUT)     :: num)
!!
!! DESCRIPTION
!!
!! Gets the number of elements in a linked list list implemented under
!! the hood.  This subroutine is overloaded and implements
!! nameValueLL_getNumReal, nameValueLL_getNumInt
!! nameValueLL_getNumStr, nameValueLL_getNumLog
!!
!! ARGUMENTS
!! context:    type of list (ie for parameters or scalars)
!! num:       number of elements in the list
!!
!!***

subroutine nameValueLL_getNumReal (context, num)

  use nameValueLL_data

  implicit none

  type (context_type), intent(in) :: context
  integer, intent(out) :: num

  num = context%n_real
end subroutine nameValueLL_getNumReal



subroutine nameValueLL_getNumInt(context, num)
  use nameValueLL_data

  implicit none
  type (context_type), intent(in) :: context
  integer, intent(out) :: num

  num = context%n_int
end subroutine nameValueLL_getNumInt




subroutine nameValueLL_getNumStr (context, num)

  use nameValueLL_data

  implicit none
  type (context_type), intent(in) :: context
  integer, intent(out) :: num

  num = context%n_str
end subroutine nameValueLL_getNumStr




subroutine nameValueLL_getNumLog(context, num)

  use nameValueLL_data
  implicit none

  type (context_type), intent(in) :: context
  integer, intent(out) :: num

  num = context%n_log
end subroutine nameValueLL_getNumLog




subroutine nameValueLL_getTotNum(context, num)

  use nameValueLL_data
  implicit none

  type (context_type), intent(in) :: context
  integer, intent(out) :: num
  integer :: count
  
  num = 0
  call  nameValueLL_getNumReal(context, count)
  num = num + count
  call nameValueLL_getNumInt(context, count)
  num = num + count
  call nameValueLL_getNumStr(context, count)
  num = num + count
  call nameValueLL_getNumLog(context, count)
  num = num + count
  return
end subroutine nameValueLL_getTotNum
