module PersonList
  use PersonNode
  use Driver, ONLY : Driver_abortFlash
  implicit none

  type person_list
     type(person_node), pointer :: H, T
  end type person_list
  
contains

  !An include file with the basic list operations AND 
  !a higher order function which operates on the nodes 
  !of the list.  The higher order function requires the 
  !following C preprocessor definition.
#define CPP_NODE_NAME person_node
#define CPP_LIST_NAME person_list
#define CPP_NODE_DEFINITION \
  use PersonNode, ONLY : person_node
#include "ut_listMethods.includeF90"

#undef CPP_NODE_NAME
#undef CPP_LIST_NAME
#undef CPP_NODE_DEFINITION

end module PersonList
