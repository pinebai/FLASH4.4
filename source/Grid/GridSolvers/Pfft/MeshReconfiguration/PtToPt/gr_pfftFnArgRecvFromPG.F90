!!****if* source/Grid/GridSolvers/Pfft/MeshReconfiguration/PtToPt/gr_pfftFnArgRecvFromPG
!!
!! NAME 
!!
!! gr_pfftFnArgRecvFromPG
!!
!! SYNOPSIS
!!
!! gr_pfftFnArgRecvFromPG(type(node), pointer  :: item)
!!
!! DESCRIPTION 
!!
!! Receives a fragment from a node in the pencil list to a node in 
!! the FLASH list.
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
!! The matching send is posted in gr_pfftFnArgSendToFG.
!! A return value of 0 indicates no error with the non-blocking send.
!!
!!***

integer function gr_pfftFnArgRecvFromPG(item)
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
       item % metadata % pfftProcID, item % metadata % tagID, & 
       comm, item % request, ierr)
  call Driver_checkMPIErrorCode(ierr)
  gr_pfftFnArgRecvFromPG = 0
end function gr_pfftFnArgRecvFromPG
