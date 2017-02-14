!!****if* source/Grid/GridSolvers/Pfft/MeshReconfiguration/Collective/gr_pfftMapFromOutput
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
!!***

subroutine gr_pfftMapFromOutput(gridVar, pfftOutputArray)
#include "constants.h"
#include "Flash.h"
#include "Pfft.h"
  use gr_pfftData, ONLY : pfft_ndim, pfft_comm
  use gr_pfftReconfigData, ONLY : pfft_maxProcData, pfft_sendBuf, &
       pfft_recvBuf
  use gr_pfftInterface, ONLY : gr_pfftBufferTransfer, &
       gr_pfftHandleJaxisFragments, gr_pfftHandleKaxisFragments
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none
  include "Flash_mpi.h"
  integer, intent(IN) :: gridVar
  real, dimension(:), target :: pfftOutputArray
  integer :: ierr, axis

  !Copy data from pfft_outArray into pfft_sendBuf.
  axis = JAXIS
  if (pfft_ndim>2) axis = KAXIS
  call gr_pfftBufferTransfer(FROM_PFFT, axis, pfft_sendBuf, pfftOutputArray)


  if (pfft_ndim>2) then

     !Perform an MPI_Alltoall to get the data in pfft_sendBuf onto 
     !the correct JAXIS processor.        
     call MPI_Alltoall(pfft_sendBuf(1,1), pfft_maxProcData, FLASH_REAL,&
          pfft_recvBuf(1,1), pfft_maxProcData, FLASH_REAL, &
          pfft_comm(IAXIS), ierr)

     !Extract data from pfft_recvBuf and place into pfft_sendBuf.
     !Data in pfft_sendBuf will be end up on the correct FLASH
     !processor after another MPI_Alltoall.
     call gr_pfftHandleKaxisFragments(FROM_PFFT)

  end if


  !Get the data onto the correct FLASH processor.
  call MPI_Alltoall(pfft_sendBuf(1,1), pfft_maxProcData, FLASH_REAL,&
       pfft_recvBuf(1,1), pfft_maxProcData, FLASH_REAL, &
       pfft_comm(IAXIS), ierr)

  !Copy data from pfft_recvBuf into the grid.
  call gr_pfftHandleJaxisFragments(FROM_PFFT, gridVar)

end subroutine gr_pfftMapFromOutput
