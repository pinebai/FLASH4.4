!!****if* source/Grid/GridSolvers/Pfft/MeshReconfiguration/PtToPt/gr_pfftMapFromOutput
!!
!! NAME 
!!
!! gr_pfftMapFromOutput
!!
!! SYNOPSIS
!!
!! gr_pfftMapFromOutput(integer (IN) :: gridVar, &
!!                      real, target :: pfftOutputArray(:))
!!
!! DESCRIPTION 
!!
!! Invokes lower level routines which copy data from the pencil output array 
!! and store it in the grid at variable "gridVar".  This is non-trivial 
!! because the data is must be placed in overlapping FLASH grid elements which 
!! may exist on different processors.
!!
!! ARGUMENTS
!!
!! gridVar - An index corresponding to the source grid variable.
!! pfftOutputArray - The output pencil array.
!!
!! SIDE EFFECTS
!!
!! Many!
!! Sets a module level pointer to pfftOutputArray.  This pointer is used 
!! by the lower level routines (which are passed as function arguments 
!! to apply_fn_to_nodes) to indirectly copy data out and into the grid.
!!
!! Sets pfft_solnGridVar to be equal to gridVar.
!!
!! Finally, the lower level routines copy data between the nodes in the 
!! FLASH fragment list (pfft_listFG) and the PFFT fragment list (pfft_listPG).
!!
!! NOTES
!!
!! 1: Copy data from the pencil array and store in nodes of the PFFT list.
!! 2: Post receives for block fragments from the nodes of the FLASH list.
!! 3: Post sends containing block fragments from the nodes of the PFFT list.
!! 4: Wait for successful delivery of messages into local FLASH list.
!! 5: Wait for successful delivery of messages into **remote** FLASH lists.
!! 6: Transfer data from nodes of the FLASH list to the grid.
!!
!! The communication between nodes is very simple, and we may actually want 
!! to slow the message exchange. This may be required if the number of messages 
!! is very large, i.e. a huge number of blocks and few processors. BUT...
!! "If the call causes some system resource to be exhausted, then it will
!! fail and return an error code.  Quality implementations of MPI should 
!! ensure that this only happens in "pathological" cases.  That is, an MPI
!! implementation should be able to support a large number of pending 
!! nonblocking operations" [MPI standard].
!! 
!!***

subroutine gr_pfftMapFromOutput(gridVar, pfftOutputArray)
  use gr_pfftReconfigData, ONLY : pfft_pfftBuf, pfft_solnGridVar, &
       pfft_listFG, pfft_listPG
  use gr_pfftNodeFnPrototypes, ONLY : gr_pfftFnArgCopyPencilToBuf, &
       gr_pfftFnArgRecvFromPG, gr_pfftFnArgSendToFG, gr_pfftFnArgMsgDelivered, &
       gr_pfftFnArgCopyBufToGrid
  use gr_pfftListObject, ONLY : apply_fn_to_nodes
  implicit none
  integer, intent(IN) :: gridVar
  real, dimension(:), target :: pfftOutputArray

  pfft_pfftBuf => pfftOutputArray
  pfft_solnGridVar = gridVar
  call apply_fn_to_nodes(gr_pfftFnArgCopyPencilToBuf, pfft_listPG) !1
  call apply_fn_to_nodes(gr_pfftFnArgRecvFromPG, pfft_listFG)      !2
  call apply_fn_to_nodes(gr_pfftFnArgSendToFG, pfft_listPG)        !3
  call apply_fn_to_nodes(gr_pfftFnArgMsgDelivered, pfft_listFG)    !4
  call apply_fn_to_nodes(gr_pfftFnArgMsgDelivered, pfft_listPG)    !5
  call apply_fn_to_nodes(gr_pfftFnArgCopyBufToGrid, pfft_listFG)   !6      
  nullify(pfft_pfftBuf)
end subroutine gr_pfftMapFromOutput
