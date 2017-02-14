module GroupList
  use GroupNode
  use Driver, ONLY : Driver_abortFlash
  implicit none

  type group_list
     type(group_node), pointer :: H, T
  end type group_list

contains

  !An include file with the basic list operations AND 
  !a higher order function which operates on the nodes 
  !of the list.  The higher order function requires the 
  !following C preprocessor definition.
#define CPP_NODE_NAME group_node
#define CPP_LIST_NAME group_list
#define CPP_NODE_DEFINITION \
  use GroupNode, ONLY : group_node
#include "ut_listMethods.includeF90"

#undef CPP_NODE_NAME
#undef CPP_LIST_NAME
#undef CPP_NODE_DEFINITION

end module GroupList
