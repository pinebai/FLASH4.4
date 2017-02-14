!!****if* source/Grid/GridSolvers/Pfft/MeshReconfiguration/PtToPt/gr_pfftMapToInput
!!
!! NAME 
!!
!! gr_pfftMapToInput
!!
!! SYNOPSIS
!!
!! gr_pfftMapToInput(integer (IN) :: gridVar, &
!!                   real, target :: pfftInputArray(:))
!!
!! DESCRIPTION 
!!
!! Invokes lower level routines which copy data from the grid at variable 
!! "gridVar" and store it into the pencil input array.  This is non-trivial 
!! because the data is obtained from overlapping FLASH grid elements which 
!! may exist on different processors.
!!
!! ARGUMENTS
!!
!! gridVar - An index corresponding to the source grid variable.
!! pfftInputArray - The input pencil array.
!!
!! SIDE EFFECTS
!!
!! Many!
!! Sets a module level pointer to pfftInputArray.  This pointer is used 
!! by the lower level routines (which are passed as function arguments 
!! to apply_fn_to_nodes) to indirectly copy data into pfftInputArray.
!!
!! Sets pfft_srcGridVar to be equal to gridVar.
!!
!! Finally, the lower level routines copy data between the nodes in the 
!! FLASH fragment list (pfft_listFG) and the PFFT fragment list (pfft_listPG).
!!
!! NOTES
!!
!! Steps:
!! 1: Copy data from the grid and store it in nodes of the FLASH list.
!! 2: Post receives for block fragments from the nodes of the PFFT list.
!! 3: Post sends containing block fragments from the nodes of the FLASH list.
!! 4: Wait for successful delivery of messages into local PFFT list.
!! 5: Wait for successful delivery of messages into **remote** PFFT lists.
!! 6: Transfer data from nodes of the PFFT list to the pencil array.
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

subroutine gr_pfftMapToInput(gridVar, pfftInputArray)
  use gr_pfftReconfigData, ONLY : pfft_pfftBuf, pfft_srcGridVar, &
       pfft_listFG, pfft_listPG
  use gr_pfftNodeFnPrototypes, ONLY : gr_pfftFnArgCopyGridToBuf, &
       gr_pfftFnArgRecvFromFG, gr_pfftFnArgSendToPG, gr_pfftFnArgMsgDelivered, &
       gr_pfftFnArgCopyBufToPencil
  use gr_pfftListObject, ONLY : apply_fn_to_nodes
  implicit none
  integer, intent(IN) :: gridVar
  real, dimension(:), target :: pfftInputArray

  pfft_pfftBuf => pfftInputArray
  pfft_srcGridVar = gridVar     
  call apply_fn_to_nodes(gr_pfftFnArgCopyGridToBuf, pfft_listFG)   !1
  call apply_fn_to_nodes(gr_pfftFnArgRecvFromFG, pfft_listPG)      !2
  call apply_fn_to_nodes(gr_pfftFnArgSendToPG, pfft_listFG)        !3
  call apply_fn_to_nodes(gr_pfftFnArgMsgDelivered, pfft_listPG)    !4
  call apply_fn_to_nodes(gr_pfftFnArgMsgDelivered, pfft_listFG)    !5
  call apply_fn_to_nodes(gr_pfftFnArgCopyBufToPencil, pfft_listPG) !6
  nullify(pfft_pfftBuf)
end subroutine gr_pfftMapToInput
