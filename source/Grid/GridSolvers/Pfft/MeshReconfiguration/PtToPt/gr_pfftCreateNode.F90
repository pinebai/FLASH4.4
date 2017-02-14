!!****if* source/Grid/GridSolvers/Pfft/MeshReconfiguration/PtToPt/gr_pfftCreateNode
!!
!! NAME 
!!
!! gr_pfftCreateNode
!!
!! SYNOPSIS
!!
!! gr_pfftCreateNode(integer (IN) :: nodeType, &
!!                   type(node), pointer  :: item)
!!
!! DESCRIPTION 
!!
!! Creates a new node on either the PFFT list or FLASH list.
!!
!! ARGUMENTS
!!
!! nodeType - 2 possibilites: PFFT_FLASH_NODE, PFFT_PFFT_NODE.
!! item - pointer to a node.  The current implementation stores the nodes in 
!!        a linked list which includes ut_listMethods.includeF90 methods.
!!
!! SIDE EFFECTS 
!!
!! Adding items to the end of the module level lists, and updating 
!! the module level counts of the number of items in each list.
!!
!! NOTES
!! 
!! We maintain two lists: one list is used to hold FLASH block fragments, 
!! and the other list is used to hold PFFT block fragments.  The lists are 
!! retained in memory so that the same communication pattern (stored in the
!! lists) can be re-used for multiple parallel FFTs.
!!
!!***

subroutine gr_pfftCreateNode(nodeType, item)
#include "Pfft.h"
  use Driver_interface, ONLY : Driver_abortFlash
  use gr_pfftNodeObject, ONLY : node, create_node
  use gr_pfftListObject, ONLY : push_back
  use gr_pfftReconfigData, ONLY : pfft_listFG, pfft_listPG, &
       pfft_numFlashNodes, pfft_numPfftNodes
  implicit none
  integer, intent(IN) :: nodeType
  type(node), pointer :: item
  
  call create_node(item)
  
  if (nodeType == PFFT_FLASH_NODE) then
     !Add the node to the end of the FLASH list.
     call push_back(pfft_listFG, item)
     pfft_numFlashNodes = pfft_numFlashNodes + 1
     
  else if (nodeType == PFFT_PENCIL_NODE) then
     !Add the node to the end of the PFFT list.
     call push_back(pfft_listPG, item)
     pfft_numPfftNodes = pfft_numPfftNodes + 1
     
  else
     call Driver_abortFlash("Node type must be send or recv!")
  end if
  
end subroutine gr_pfftCreateNode
