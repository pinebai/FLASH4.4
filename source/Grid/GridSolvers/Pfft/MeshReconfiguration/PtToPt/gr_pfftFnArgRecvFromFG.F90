!!****if* source/Grid/GridSolvers/Pfft/MeshReconfiguration/PtToPt/gr_pfftFnArgRecvFromFG
!!
!! NAME 
!!
!! gr_pfftFnArgRecvFromFG
!!
!! SYNOPSIS
!!
!! gr_pfftFnArgRecvFromFG(type(node), pointer  :: item)
!!
!! DESCRIPTION 
!!
!! Receives a fragment from a node in the FLASH list to a node in 
!! the pencil node. 
!!
!! ARGUMENTS
!!
!! item - pointer to a node.  The current implementation stores the nodes in 
!!        a linked list which includes ut_listMethods.includeF90 methods.
!!
!! SIDE EFFECTS 
!!
!! A non-blocking MPI communication.
!!
!! NOTES
!! 
!! The matching send is posted in gr_pfftFnArgSendToPG.
!! A return value of 0 indicates no error with the non-blocking receive.
!!
!!***

integer function gr_pfftFnArgRecvFromFG(item)
#include "constants.h"
  use Driver_interface, ONLY : Driver_checkMPIErrorCode
  use gr_pfftNodeObject, ONLY : node
  use gr_pfftData, ONLY : pfft_comm

  implicit none
  type(node), pointer  :: item
  include "Flash_mpi.h"
  integer :: ierr, comm

  comm = pfft_comm(IAXIS)
  call MPI_Irecv(item % buf(1), size(item % buf(:)), FLASH_REAL, &
       item % metadata % flashProcID, item % metadata % tagID, & 
       comm, item % request, ierr)
  call Driver_checkMPIErrorCode(ierr)
  gr_pfftFnArgRecvFromFG = 0
end function gr_pfftFnArgRecvFromFG
