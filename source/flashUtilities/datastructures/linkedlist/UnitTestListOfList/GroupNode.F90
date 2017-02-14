module GroupNode
  use PersonList, ONLY : person_initialise_list => initialise_list, &
       person_finalise_list => finalise_list, &
       person_print_list => print_list, person_list
  use Driver, ONLY : Driver_abortFlash
  implicit none

  type group_node
     character (len=80) :: category
     type(person_list), pointer :: person_list
     type(group_node), pointer :: next, prev
  end type group_node

contains

  !The user must provide: create_node, destroy_node and print_node.
  subroutine create_node(item)
    implicit none
    type(group_node), pointer :: item
    integer :: err

    allocate(item, STAT=err)
    if (err /= 0) then
       call Driver_abortFlash ("[GroupNode::create_node]: "//&
            "Memory cannot be allocated")
    end if
    item % category = "NULL"
    call person_initialise_list(item % person_list)
    nullify(item % next, item % prev)
  end subroutine create_node

  subroutine destroy_node(item)
    implicit none
    type(group_node), pointer  :: item
    integer :: err

    call person_finalise_list(item % person_list)
    nullify(item % person_list)

    nullify(item % next, item % prev)
    deallocate(item, STAT=err)
    if (err /= 0) then
       call Driver_abortFlash ("[GroupNode::destroy_node]: "//&
            "Memory cannot be deallocated")
    end if
    nullify(item)
  end subroutine destroy_node

  subroutine print_node(unitNumber, item)    
    implicit none
    integer, intent(IN) :: unitNumber
    type(group_node), pointer :: item

    write(unitNumber,'(2A)') "Category: ", trim(item % category)
    call person_print_list(item % person_list, unitNumber)
    write(unitNumber,*) ""
  end subroutine print_node
end module GroupNode
