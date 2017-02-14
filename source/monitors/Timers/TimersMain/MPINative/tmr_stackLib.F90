!!****if* source/monitors/Timers/TimersMain/MPINative/tmr_stackLib
!!
!! NAME
!!
!!   tmr_stackLib
!!
!! SYNOPSIS
!!   tmr_stackLib(tmr_stack(INOUT) :: instack)
!!
!! DESCRIPTION
!!   This is a collection of routines to work on stacks of integers
!!   (the routines prefaced by tmr_stack) and routines to work on lists
!!   of those stacks (the routines prefaced by tmr_stackList).  The
!!   routines generally take arguments of type tmr_stack or
!!   tmr_stackList, and these types are defined in Timers_data.
!!
!! ARGUMENTS
!!
!! instack :
!!
!!***

subroutine tmr_stackZero(instack)

  use Timers_data, ONLY: tmr_stack

  implicit none

  type (tmr_stack), intent(INOUT) :: instack

  instack%stackPointer = 0
  instack%stack = 0

end subroutine tmr_stackZero

subroutine tmr_stackInit(instack, list, size, result)

  use Timers_data, ONLY: tmr_stack

  implicit none

  type(tmr_stack), intent(OUT) :: instack
  integer, intent(IN) :: size
  integer, intent(IN) :: list(size)
  integer, intent(out) :: result

  integer :: i
  
  call tmr_stackZero(instack)
  result = 1
  do i = 1, size
     call tmr_stackPush(instack, list(i), result)
     if (result < 0) then
        exit
     endif
  end do

end subroutine tmr_stackInit

subroutine tmr_stackLength(instack, length)

  use Timers_data, ONLY: tmr_stack

  implicit none

  type(tmr_stack), intent(IN) :: instack
  integer, intent(out) :: length
  
  length = instack%stackPointer

  return

end subroutine tmr_stackLength

subroutine tmr_stackAssign(receiver, giver)
  use Timers_data, ONLY: tmr_Stack

  implicit none

  type(tmr_stack), intent(INOUT) :: receiver
  type(tmr_stack), intent(IN)    :: giver
  
  receiver%stack = giver%stack
  receiver%stackPointer = giver%stackPointer

end subroutine tmr_stackAssign

subroutine tmr_stackPush(instack, i, ierr)

  use Timers_data, ONLY: tmr_stack, tmr_maxCallStackDepth

  implicit none

  type(tmr_stack), intent(INOUT):: instack
  integer, intent(IN) :: i
  integer, intent(out) :: ierr !this is the error type
  
  if (instack%stackPointer < tmr_maxCallStackDepth) then
     instack%stackPointer = instack%stackPointer + 1
     instack%stack(instack%stackPointer) = i
     ierr = instack%stackPointer
  else
     ierr = -1
  end if

end subroutine tmr_stackPush

subroutine tmr_stackPop(instack, ierr)

  use Timers_data, ONLY: tmr_stack

  implicit none

  integer, intent(out)         :: ierr
  type(tmr_stack), intent(inout)  :: instack
  
  if (instack%stackPointer > 0) then
     ierr = instack%stack(instack%stackPointer)
     instack%stackPointer = instack%stackPointer - 1
  else
     ierr = -1
  end if

end subroutine tmr_stackPop

subroutine tmr_stackTop(instack, val)
 
  use Timers_data, ONLY: tmr_stack

  implicit none

  integer, intent(out)        :: val
  type(tmr_stack), intent(in) :: instack
  
  if (instack%stackPointer > 0) then
     val = instack%stack(instack%stackPointer)
  else
     val = -1
  end if
  return

end subroutine tmr_stackTop

subroutine tmr_stacksEqual(stack1, stack2, bool)

  use Timers_data, ONLY: tmr_stack

  implicit none

  logical, intent(out) :: bool
  type (tmr_stack), intent(IN) :: stack1, stack2
  integer i
  
  bool = .true.
  if (stack1%stackPointer == stack2%stackPointer) then
     do i = 1, stack1%stackPointer
        if (stack1%stack(i) .ne. stack2%stack(i)) then
           bool = .false.
           exit
        end if
     end do
  else
     bool = .false.
  end if

  return

end subroutine tmr_stacksEqual

subroutine tmr_stackGet (instack, index, val)

  use Timers_data, ONLY: tmr_stack

  implicit none
  type(tmr_stack), intent(in)  :: instack
  integer, intent(in)          :: index
  integer, intent(out)         :: val
  
  val = instack%stack(index)
  return 

end subroutine tmr_stackGet

subroutine tmr_stackListZero(instackList)

  use Timers_data, ONLY: tmr_stackList

  implicit none

  type(tmr_stackList), intent(inout)  :: instackList
  
  instackList%lastItem = 0

end subroutine tmr_stackListZero

subroutine tmr_stackListsEqual(instackList1, instackList2, bool)

  use Timers_data, ONLY: tmr_stackList

  implicit none

  type(tmr_stackList), intent(in)  :: instackList1, instackList2
  logical, intent(out) :: bool

  logical stacksSame
  integer i

  bool = .true.
  stacksSame = .false.

  if (instackList1%lastItem .eq. instackList2%lastItem) then
     do i = 1, instackList1%lastItem
        call tmr_stacksEqual(instackList1%stacks(i), instackList2%stacks(i), stacksSame)
        if (.not. stacksSame) then
           bool = .false.
           exit
        end if
     end do
  else
     bool = .false.
  end if

  return
end subroutine tmr_stackListsEqual

subroutine tmr_stackListAdd(inlist, instack, val)

  use Timers_data, ONLY: tmr_stack, tmr_stackList, tmr_maxTimerParents

  implicit none

  integer, intent(out) :: val
  type(tmr_stackList), intent(INOUT) :: inlist
  type(tmr_stack),  intent(IN)        :: instack
  
  if (inlist%lastItem == tmr_maxTimerParents) then
     val = -1
  else 
     inlist%lastItem = inlist%lastItem + 1
     call tmr_stackAssign(inlist%stacks(inlist%lastItem), instack)
     val = inlist%lastItem
  end if

  return

end subroutine tmr_stackListAdd

subroutine tmr_stackListGet(instackList, instack, index)

  use Timers_data, ONLY: tmr_stack, tmr_stackList

  implicit none

  type(tmr_stackList), intent(IN) :: instackList
  type(tmr_stack), intent(INOUT) :: instack
  integer, intent(IN) :: index
  
  call tmr_stackAssign(instack, instackList%stacks(index))

end subroutine tmr_stackListGet

subroutine tmr_stackListLength(instackList, length)

  use Timers_data, ONLY: tmr_stackList

  implicit none

  type(tmr_stackList), intent(IN) :: instackList
  integer, intent(out) ::  length
  
  length = instackList%lastItem

  return

end subroutine tmr_stackListLength


subroutine tmr_stackListIndex(instackList, instack, ret)

  use Timers_data, ONLY: tmr_stack, tmr_stackList

  implicit none

  type(tmr_stackList), intent(IN) :: instackList
  type(tmr_stack), intent(IN) :: instack
  integer, intent(out) :: ret

  logical :: bool

  integer i
  ret = 0
  if (instackList%lastItem > 0) then
     do i = 1, instackList%lastItem
        call tmr_stacksEqual(instackList%stacks(i),instack, bool)
        if (bool) then
           ret = i 
           exit
        end if
     end do
  end if

  return
end subroutine tmr_stackListIndex

