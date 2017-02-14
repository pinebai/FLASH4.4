Program UnitTest

  use UnitTest_ListObject
  implicit none

  interface
     integer function fnPrintMagicValue(item)
       use UnitTest_NodeObject, only : node
       implicit none
       type(node), pointer  :: item
     end function fnPrintMagicValue
  end interface
  interface
     logical function fnFindMagicValue(item)
       use UnitTest_NodeObject, only : node  
       implicit none
       type(node), pointer  :: item
     end function fnFindMagicValue
  end interface
  interface
     integer function fnLongFunction(item)
       use UnitTest_NodeObject, only : node  
       implicit none
       type(node), pointer  :: item
     end function fnLongFunction
  end interface
  interface
     subroutine CreateNodeWithValue(item, value)
       use UnitTest_NodeObject, only : node, create_node
       implicit none
       type(node), pointer :: item
       integer, intent(IN) :: value
     end subroutine CreateNodeWithValue
  end interface
  interface
     subroutine CheckPointerState(L)
       use UnitTest_ListObject, only : list
       implicit none
       type(list), pointer :: L
     end subroutine CheckPointerState
  end interface


  type(list), pointer :: L
  type(node), pointer :: item
  integer, parameter :: Nitems = 4, unitStdout = 6
  integer :: i


  call initialise_list(L)


  !TEST: Push nodes to front of list and delete from front
  print *, "------------------- TEST1 ------------------------"
  do i = 1, Nitems
     call CreateNodeWithValue(item, i)
     call push_front(L, item)
  end do
  do i = 1, Nitems+1
     call print_list(L, unitStdout)
     call pop_front(L)
  end do
  call CheckPointerState(L)


  !TEST: Push nodes to front of list and delete from back
  print *, ""
  print *, ""
  print *, "------------------- TEST2 ------------------------"

  do i = 1, Nitems
     call CreateNodeWithValue(item, i)
     call push_front(L, item)
  end do
  do i = 1, Nitems+1
     call print_list(L, unitStdout)
     call pop_back(L)
  end do
  call CheckPointerState(L)


  !TEST: Push nodes to back of list and delete from back
  print *, ""
  print *, ""
  print *, "------------------- TEST3 ------------------------"
  do i = 1, Nitems
     call CreateNodeWithValue(item, i)
     call push_back(L, item)
  end do
  do i = 1, Nitems+1
     call print_list(L, unitStdout)
     call pop_back(L)
  end do
  call CheckPointerState(L)


  !TEST: Push nodes to back of list and delete from front
  print *, ""
  print *, ""
  print *, "------------------- TEST4 ------------------------"
  do i = 1, Nitems
     call CreateNodeWithValue(item, i)
     call push_back(L, item)
  end do
  do i = 1, Nitems+1
     call print_list(L, unitStdout)
     call pop_front(L)
  end do
  call CheckPointerState(L)



  !TEST: Do a bit of a random mix.
  print *, ""
  print *, ""
  print *, "------------------- TEST5 ------------------------"
  do i = 1, 3
     call CreateNodeWithValue(item, i)
     call push_back(L, item)
  end do
  call pop_back(L)
  call pop_back(L)
  call print_list(L, unitStdout) !1

  do i = 4, 8
     call CreateNodeWithValue(item, i)
     call push_front(L, item)
  end do
  call print_list(L, unitStdout) !8,7,6,5,4,1

  call pop_front(L)
  call pop_front(L)
  call pop_front(L)
  call print_list(L, unitStdout) !5,4,1

  call pop_back(L)
  call pop_back(L)
  call print_list(L, unitStdout) !5


  call pop_front(L)
  call print_list(L, unitStdout) !empty

  call pop_front(L)
  call print_list(L, unitStdout) !empty
  call CheckPointerState(L)



  !TEST: Delete list in one go (back to front): Deallocates list container too.
  print *, ""
  print *, ""
  print *, "------------------- TEST6 ------------------------"
  do i = 1, 6
     call CreateNodeWithValue(item, i)
     call push_back(L, item)
  end do
  call print_list(L, unitStdout)
  call finalise_list(L, fromFront=.false.)
  call print_list(L, unitStdout)




  !TEST: Delete list in one go (front to back): Deallocates list container too.
  print *, ""
  print *, ""
  print *, "------------------- TEST7 ------------------------"
  call initialise_list(L) !We earlier called finalise, so must re-initialise.
  do i = 1, 6
     call CreateNodeWithValue(item, i)
     call push_back(L, item)
  end do
  call print_list(L, unitStdout)
  call finalise_list(L, fromFront=.true.)
  call print_list(L, unitStdout)



  !Tests the higher order functions: apply_fn_to_nodes & get_matching_node
  print *, ""
  print *, ""
  print *, "------------------- TEST8 ------------------------"
  call initialise_list(L) !We earlier called finalise, so must re-initialise.
  do i = 1, 4
     call CreateNodeWithValue(item, i)
     call push_back(L, item)
  end do

  call apply_fn_to_nodes(fnPrintMagicValue, L)

  !Search in direction HEAD to TAIL for matching node.
  call get_matching_node(fnFindMagicValue, L % H, item, .true.)
  print *, ""
  print *, ""
  if (.not.associated(item)) then
     print *, "ERROR matching node not found"
  else
     print *, "H->T This is the matching node from our search:"
     call print_node(unitStdout, item)
  end if

  !Search in direction TAIL to HEAD for matching node.
  call get_matching_node(fnFindMagicValue, L % T, item, .false.)
  print *, ""
  print *, ""
  if (.not.associated(item)) then
     print *, "ERROR matching node not found"
  else
     print *, "T->H This is the matching node from our search:"
     call print_node(unitStdout, item)
  end if
  print *, ""
  print *, ""
  call print_list(L, unitStdout)


  print *, ""
  print *, "HERE COMES THE THREADED SECTION"
  call map_impure_fn(fnLongFunction, L, headToTailDirection=.true., parallelFlag=.true.)
  call print_list(L, unitStdout)
  print *, "FINISHED THREADED SECTION"
  print *, ""

  call finalise_list(L)  !Same as finalise_list(L, fromFront=.false.)
  call print_list(L, unitStdout)

End Program UnitTest


subroutine CreateNodeWithValue(item, value)
  use UnitTest_NodeObject, only : node, create_node
  implicit none
  type(node), pointer :: item
  integer, intent(IN) :: value

  call create_node(item)
  item % SomeInt = value
end subroutine CreateNodeWithValue


subroutine CheckPointerState(L)
  use UnitTest_ListObject, only : list
  implicit none
  type(list), pointer :: L

  if (associated(L % H)) then
     print *, "Head pointer not null!"
     stop
  end if
  if (associated(L % T)) then
     print *, "Tail pointer not null!"
     stop
  end if
end subroutine CheckPointerState


integer function fnPrintMagicValue(item)
  use UnitTest_NodeObject, only : node  
  implicit none
  type(node), pointer  :: item
  integer, parameter :: MagicValue = 2
  if (item % SomeInt == MagicValue) then
     print *, ":)        Found the magic value:", MagicValue
  else
     print *, ":( Did not find the magic value:", MagicValue
  end if
  fnPrintMagicValue = 0
end function fnPrintMagicValue


logical function fnFindMagicValue(item)
  use UnitTest_NodeObject, only : node  
  implicit none
  type(node), pointer  :: item
  integer, parameter :: MagicValue = 2
  if (item % SomeInt == MagicValue) then
     fnFindMagicValue = .true.
  else
     fnFindMagicValue = .false.
  end if
end function fnFindMagicValue


integer function fnLongFunction(item)
  !$ use omp_lib
  use UnitTest_NodeObject, only : node  
  implicit none
  type(node), pointer  :: item
  integer :: threadThatTouched

  threadThatTouched = -1
  !$ threadThatTouched = omp_get_thread_num()

  item % threadThatTouched = threadThatTouched
  call sleep(1)
  fnLongFunction = 0
end function fnLongFunction
