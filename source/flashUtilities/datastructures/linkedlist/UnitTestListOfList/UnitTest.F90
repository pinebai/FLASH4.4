!This is just a small test program that demonstrates how to create 
!a generic list of lists.  Pfft unit uses a list of lists, so it 
!is nice to test the list code independent of FLASH.

Program UnitTest

  use PersonList, ONLY : person_list, person_initialise_list => initialise_list, &
       person_finalise_list => finalise_list, person_print_list => print_list, &
       person_push_back => push_back, person_pop_front => pop_front, person_node, &
       person_create_node => create_node
  use GroupList, ONLY : group_list, group_initialise_list => initialise_list, &
       group_finalise_list => finalise_list, group_print_list => print_list, &
       group_push_back => push_back, group_pop_front => pop_front, group_node, &
       group_create_node => create_node
  implicit none
  type(person_node), pointer :: person
  type(group_list), pointer :: Lgroup
  type(group_node), pointer :: group
  integer, parameter :: unitStdout = 6


  call group_initialise_list(Lgroup)


  !Initialise goalkeepers
  call group_create_node(group)
  group % category = "Goalkeeper"
  call group_push_back(Lgroup, group)

  call person_create_node(person)
  person % name = "Jon Busch"
  person % age = 27
  call person_push_back(group % person_list, person)


  !Initialise defenders
  call group_create_node(group)
  group % category = "Defenders"
  call group_push_back(Lgroup, group)

  call person_create_node(person)
  person % name = "Tim Ward"
  person % age = 25
  call person_push_back(group % person_list, person)

  call person_create_node(person)
  person % name = "John Thorrington"
  person % age = 28
  call person_push_back(group % person_list, person)


  !Initialise midfielders
  call group_create_node(group)
  group % category = "Midfielders"
  call group_push_back(Lgroup, group)

  call person_create_node(person)
  person % name = "Marco Pappa"
  person % age = 23
  call person_push_back(group % person_list, person)

  call person_create_node(person)
  person % name = "Logan Pause"
  person % age = 31
  call person_push_back(group % person_list, person)


  !Print List of lists.
  call group_print_list(Lgroup, unitStdout)

  !Free all memory.
  call group_finalise_list(Lgroup)

End Program UnitTest
