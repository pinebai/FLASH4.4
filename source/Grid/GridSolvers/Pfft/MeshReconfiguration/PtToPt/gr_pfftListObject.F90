!!****if* source/Grid/GridSolvers/Pfft/MeshReconfiguration/PtToPt/gr_pfftListObject
!!
!! NAME
!!  gr_pfftListObject
!!
!! SYNOPSIS
!!
!!  use gr_pfftListObject
!!
!! 
!!   
!!***
module gr_pfftListObject
  use gr_pfftNodeObject
  implicit none

  type list
     type(node), pointer :: H, T
  end type list

contains

  !An include file with the list operations.
#define CPP_NODE_DEFINITION \
  use gr_pfftNodeObject, ONLY : node
#include "ut_listMethods.includeF90"
#undef CPP_NODE_DEFINITION
end module gr_pfftListObject
