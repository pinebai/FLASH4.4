!!****ih* source/Grid/GridSolvers/Pfft/MeshReconfiguration/PtToPt/gr_pfftNodeFnPrototypes
!!
!! NAME 
!!   gr_pfftNodeFnPrototypes
!!
!! SYNOPSIS
!!   
!!   use gr_pfftNodeFnPrototypes
!!
!! DESCRIPTION
!! 
!!  This is the header file for procedures which use the node datatype.
!!  We have placed the prototypes here because they depend on an 
!!  external datatype which is only created in PFFT point to point sub-unit.
!!
!!***

module gr_pfftNodeFnPrototypes
  implicit none

  interface
     integer function gr_pfftFnArgRecvFromFG(item)
       use gr_pfftNodeObject, ONLY : node
       implicit none
       type(node), pointer  :: item
     end function gr_pfftFnArgRecvFromFG
  end interface

  interface
     integer function gr_pfftFnArgSendToPG(item)
       use gr_pfftNodeObject, ONLY : node
       implicit none
       type(node), pointer  :: item
     end function gr_pfftFnArgSendToPG
  end interface

  interface
     integer function gr_pfftFnArgRecvFromPG(item)
       use gr_pfftNodeObject, ONLY : node
       implicit none
       type(node), pointer  :: item
     end function gr_pfftFnArgRecvFromPG
  end interface

  interface
     integer function gr_pfftFnArgSendToFG(item)
       use gr_pfftNodeObject, ONLY : node
       implicit none
       type(node), pointer  :: item
     end function gr_pfftFnArgSendToFG
  end interface

  interface
     integer function gr_pfftFnArgMsgDelivered(item)
       use gr_pfftNodeObject, ONLY : node
       implicit none
       type(node), pointer  :: item
     end function gr_pfftFnArgMsgDelivered
  end interface

  interface
     integer function gr_pfftFnArgCopyGridToBuf(item)
       use gr_pfftNodeObject, ONLY : node
       implicit none
       type(node), pointer  :: item
     end function gr_pfftFnArgCopyGridToBuf
  end interface

  interface
     integer function gr_pfftFnArgCopyBufToGrid(item)
       use gr_pfftNodeObject, ONLY : node
       implicit none
       type(node), pointer  :: item
     end function gr_pfftFnArgCopyBufToGrid
  end interface

  interface
     integer function gr_pfftFnArgCopyPencilToBuf(item)
       use gr_pfftNodeObject, ONLY : node
       implicit none
       type(node), pointer  :: item
     end function gr_pfftFnArgCopyPencilToBuf
  end interface

  interface
     integer function gr_pfftFnArgCopyBufToPencil(item)
       use gr_pfftNodeObject, ONLY : node
       implicit none
       type(node), pointer  :: item
     end function gr_pfftFnArgCopyBufToPencil
  end interface

  interface
     subroutine gr_pfftCreateNode(nodeType, item)
       use gr_pfftNodeObject, ONLY : node
       implicit none
       integer, intent(IN) :: nodeType
       type(node), pointer :: item
     end subroutine gr_pfftCreateNode
  end interface
end module gr_pfftNodeFnPrototypes
