!!****if* source/Grid/GridSolvers/Pfft/MeshReconfiguration/Collective/gr_pfftMapToInput
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
!!***

subroutine gr_pfftMapToInput(gridVar, pfftInputArray)
#include "constants.h"
#include "Flash.h"
#include "Pfft.h"
  use gr_pfftData, ONLY : pfft_ndim, pfft_comm
  use gr_pfftReconfigData, ONLY : pfft_sendBuf, pfft_recvBuf, pfft_maxProcData
  use gr_pfftInterface, ONLY : gr_pfftBufferTransfer, &
       gr_pfftHandleJaxisFragments, gr_pfftHandleKaxisFragments
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none
  integer, intent(IN) :: gridVar
  real, dimension(:), target :: pfftInputArray
  include "Flash_mpi.h"
  integer :: ierr, axis

  !This routine will copy the values from the grid (i.e. held in
  !solnData(gridVar,i,j,k)) and store them in pfft_sendBuf.  The 
  !data in pfft_sendBuf is ordered so that after an MPI_Alltoall it 
  !will end up on a "PFFT" processor.  Depending on the dimensionality 
  !of the problem a further communication between PFFT processors 
  !may be neccessary.
  call gr_pfftHandleJaxisFragments(TO_PFFT, gridVar)

  call MPI_Alltoall(pfft_sendBuf(1,1), pfft_maxProcData, FLASH_REAL,&
       pfft_recvBuf(1,1), pfft_maxProcData, FLASH_REAL, &
       pfft_comm(IAXIS), ierr)


  if (pfft_ndim > 2) then

     !In a 3D problem we need to handle KAXIS data movement too.
     !Here, we copy data from pfft_recvBuf into pfft_sendBuf, 
     !so that after the MPI_Alltoall, the data is on the 
     !correct PFFT processor.
     call gr_pfftHandleKaxisFragments(TO_PFFT)

     call MPI_Alltoall(pfft_sendBuf(1,1), pfft_maxProcData, FLASH_REAL,&
          pfft_recvBuf(1,1), pfft_maxProcData, FLASH_REAL, &
          pfft_comm(IAXIS), ierr)

  endif


  !The data in pfft_recvBuf is on the correct processor.  Therefore, the 
  !final step is to copy data from pfft_recvBuf into pfft_inArray.
  axis = JAXIS
  if (pfft_ndim > 2) axis = KAXIS
  call gr_pfftBufferTransfer(TO_PFFT, axis, pfft_recvBuf, pfftInputArray)

end subroutine gr_pfftMapToInput
