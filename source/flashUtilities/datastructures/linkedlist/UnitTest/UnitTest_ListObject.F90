module UnitTest_ListObject
  use UnitTest_NodeObject
  implicit none

  type list
     type(node), pointer :: H, T
  end type list

contains

  !An include file with the basic list operations AND 
  !a higher order function which operates on the nodes 
  !of the list.  The higher order function requires the 
  !following C preprocessor definition.
#define CPP_NODE_DEFINITION \
  use UnitTest_NodeObject
#include "ut_listMethods.includeF90"
#undef CPP_NODE_DEFINITION
end module UnitTest_ListObject
