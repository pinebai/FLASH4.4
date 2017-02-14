module PersonNode
  use Driver, ONLY : Driver_abortFlash
  implicit none

  type person_node
     character (len=80) :: name
     integer :: age
     type(person_node), pointer :: next, prev
  end type person_node

contains

  !The user must provide: create_node, destroy_node and print_node.
  subroutine create_node(item)
    implicit none
    type(person_node), pointer :: item
    integer :: err    
    allocate(item, STAT=err)
    if (err /= 0) then
       call Driver_abortFlash ("[PersonNode::create_node]: "//&
            "Memory cannot be allocated")
    end if
    item % name = "NULL"
    item % age = -1
    nullify(item % next, item % prev)
  end subroutine create_node


  subroutine destroy_node(item)
    implicit none
    type(person_node), pointer  :: item
    integer :: err

    deallocate(item, STAT=err)
    if (err /= 0) then
       call Driver_abortFlash ("[PersonNode::destroy_node]: "//&
            "Memory cannot be deallocated")
    end if
    nullify(item)
  end subroutine destroy_node


  subroutine print_node(unitNumber, item)    
    implicit none
    integer, intent(IN) :: unitNumber
    type(person_node), pointer :: item
    write(unitNumber,'(3A,i3)') "Name: ", trim(item % name), &
         ", age: ", item % age
  end subroutine print_node
end module PersonNode
