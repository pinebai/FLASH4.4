!!****if* source/Grid/GridMain/paramesh/paramesh4/gr_freeCommRecvBuffer
!!
!! NAME
!!
!!  gr_freeCommRecvBuffer
!!
!! SYNOPSIS
!!
!!  gr_freeCommRecvBuffer()
!!
!! DESCRIPTION
!!
!!  Free the communication buffer temprecv_buf that is allocated by the PARAMESH
!!  subroutine mpi_amr_comm_setup.
!!
!!  This buffer is used for receiving MPI messages received during processing
!!  of guard cells filling, restriction, prolongation, and flux / edge correction.
!!  Normally, when not in no_permanent_guardcells mode, the buffer is not needed
!!  any more when control after these operations returns to FLASH.  Since sizes
!!  of the buffer can be very large, it may be worth freeing it in between uses.
!!  When this is not done, it will be freed anyway (and then allocated again,
!!  likely with a different size) the next time mpi_amr_comm_setup is called.
!!
!!
!! ARGUMENTS
!!
!!   No arguments
!!
!! NOTES
!!
!!  In general, this subroutine should not be called by FLASH in
!!  no_permanent_guardcells mode: In that mode, the buffer needs to
!!  stay around so it can be used as a source of remote cell data
!!  when Grid_getBlkPtr is called.
!!
!! HISTORY
!!  Created gr_freeCommRecvBuffer   - KW Jan 2009
!!***


subroutine gr_freeCommRecvBuffer()
  use mpi_morton, ONLY: temprecv_buf
  implicit none

  if(allocated(temprecv_buf)) deallocate(temprecv_buf)
end subroutine gr_freeCommRecvBuffer
