module UnitTest_NodeObject

  implicit none
  type node
     integer :: someInt, threadThatTouched
     type(node), pointer :: next, prev
  end type node

  contains

  !The user must provide: create_node, destroy_node and print_node.
  subroutine create_node(item)
    implicit none
    type(node), pointer :: item
    integer :: err    
    allocate(item, STAT=err)
    if (err /= 0) then
       call Driver_abortFlash &
            ("[create_node]: No more heap memory")
    end if
    item % SomeInt = -1
    item % threadThatTouched = -1
    nullify(item % next, item % prev)
  end subroutine create_node


  subroutine destroy_node(item)
    implicit none
    type(node), pointer  :: item    
    deallocate(item); nullify(item)
  end subroutine destroy_node


  subroutine print_node(unitNumber, item)    
    implicit none
    integer, intent(IN) :: unitNumber
    type(node), pointer :: item
    write(unitNumber,*) "[print_node]: Value", item % SomeInt, &
         "(Optional Thread)", item % threadThatTouched
  end subroutine print_node


  !Just whack this here, so we can test without FLASH.
  subroutine Driver_abortFlash(msg)
    implicit none
    character (len=*), intent(IN) :: msg
    print *, "ERROR!!!"
    print *, msg
    stop
  end subroutine Driver_abortFlash

end module UnitTest_NodeObject
