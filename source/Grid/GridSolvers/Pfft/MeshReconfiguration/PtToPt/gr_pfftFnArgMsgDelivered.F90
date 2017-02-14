!!****if* source/Grid/GridSolvers/Pfft/MeshReconfiguration/PtToPt/gr_pfftFnArgMsgDelivered
!!
!! NAME 
!!
!! gr_pfftFnArgMsgDelivered
!!
!! SYNOPSIS
!!
!! gr_pfftFnArgMsgDelivered(type(node), pointer  :: item)
!!
!! DESCRIPTION 
!!
!! Waits for a message to be delivered into passed node, or for a 
!! message to be delivered into a remote node.
!!
!! ARGUMENTS
!!
!! item - pointer to a node.  The current implementation stores the nodes in 
!!        a linked list which includes ut_listMethods.includeF90 methods.
!!
!! SIDE EFFECTS 
!!
!! A blocking MPI_Wait.
!!
!! NOTES
!! 
!! A return value of 0 indicates no error with message delivery.
!!
!!***

integer function gr_pfftFnArgMsgDelivered(item)
  use Driver_interface, ONLY : Driver_checkMPIErrorCode
  use gr_pfftNodeObject, ONLY : node

  implicit none
  type(node), pointer  :: item
  include "Flash_mpi.h"
  integer :: ierr

  call MPI_Wait(item % request, item % status, ierr)
  call Driver_checkMPIErrorCode(ierr)
  gr_pfftFnArgMsgDelivered = 0
end function gr_pfftFnArgMsgDelivered
