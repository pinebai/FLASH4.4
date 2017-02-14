!!****if* source/Grid/GridSolvers/Pfft/Grid_pfftFinalize
!!
!! NAME 
!!
!!   Grid_pfftFinalize
!!
!! SYNOPSIS
!!
!!   Grid_pfftFinalize()
!!
!! DESCRIPTION 
!!
!!  Deallocate the work array and trignometric tables used in transforms
!!
!!***

subroutine Grid_pfftFinalize()
#include "constants.h"
  use gr_pfftData, ONLY : pfft_trigIaxis, pfft_trigJaxis, pfft_trigKaxis, &
       pfft_work1, pfft_work2, pfft_needMap, pfft_comm, pfft_ndim, &
       pfft_wave, pfft_usableProc, pfft_commWithTopology, &
       pfft_isSimpleNonMappedAMR1D
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_data, ONLY : gr_meshComm

  implicit none
  include "Flash_mpi.h"
  integer :: ierr, error

  if (.not.pfft_usableProc) return


  if (pfft_needMap .or. pfft_isSimpleNonMappedAMR1D) then
     !If the IAXIS communicator has been split, we must free the communicator.
     if (pfft_comm(IAXIS) /= gr_meshComm) then
        call MPI_Comm_free(pfft_comm(IAXIS), ierr)
     end if

     call MPI_Comm_free(pfft_comm(JAXIS), ierr)
     call MPI_Comm_free(pfft_comm(KAXIS), ierr)
     call MPI_Comm_free(pfft_commWithTopology, ierr)
  end if


  if (pfft_needMap) then
     call gr_pfftFreeMap()
  end if


  deallocate(pfft_trigIaxis, STAT=error)
  if (error /= 0) then
     call Driver_abortFlash("[Grid_pfftFinalize]: " // &
          "Severe error. pfft_trigIaxis cannot be deallocated!")
  end if

  deallocate(pfft_wave, STAT=error)
  if (error /= 0) then
     call Driver_abortFlash("[Grid_pfftFinalize]: " // &
          "Severe error. pfft_wave cannot be deallocated!")
  end if

  deallocate(pfft_work1, pfft_work2, STAT=error)
  if (error /= 0) then
     call Driver_abortFlash("[Grid_pfftFinalize]: " // &
          "Severe error. work arrays cannot be deallocated!")
  end if

  if(pfft_ndim > 1) then
     deallocate(pfft_trigJaxis, STAT=error)
     if (error /= 0) then
        call Driver_abortFlash("[Grid_pfftFinalize]: " // &
             "Severe error. pfft_trigJaxis cannot be deallocated!")
     end if

  end if

  if(pfft_ndim > 2) then
     deallocate(pfft_trigKaxis, STAT=error)
     if (error /= 0) then
        call Driver_abortFlash("[Grid_pfftFinalize]: " // &
             "Severe error. pfft_trigKaxis cannot be deallocated!")
     end if
  end if
end subroutine Grid_pfftFinalize
