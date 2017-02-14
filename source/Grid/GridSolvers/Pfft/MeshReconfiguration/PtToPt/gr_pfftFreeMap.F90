!!****if* source/Grid/GridSolvers/Pfft/MeshReconfiguration/PtToPt/gr_pfftFreeMap
!!
!! NAME 
!!
!! gr_pfftFreeMap
!!
!! SYNOPSIS
!!
!! gr_pfftFreeMap()
!!
!! DESCRIPTION 
!!
!! Finalises variables, arrays and data structure which were used for
!! mapping from FLASH grid -> Pencil grid and Pencil grid -> FLASH grid. 
!!
!! ARGUMENTS
!!
!! SIDE EFFECTS 
!!
!! Finalisation of module level variables, arrays, lists:
!!
!! NOTES
!! 
!! We destroy an MPI datatype in this routine which describes the metadata 
!! that we exchange between nodes.
!!
!!***

subroutine gr_pfftFreeMap()
#include "constants.h"
  use Driver_interface, ONLY : Driver_checkMPIErrorCode, &
    Driver_abortFlash
  use Logfile_interface, ONLY : Logfile_close
  use gr_pfftListObject, ONLY : finalise_list
  use gr_pfftData, ONLY : pfft_ndim
  use gr_pfftReconfigData, ONLY : pfft_metaType, pfft_numFlashNodes, &
       pfft_numPfftNodes, pfft_srcGridVar, pfft_solnGridVar, pfft_pfftBuf, &
       pfft_numMsgSendToEachProc, pfft_listFG, pfft_listPG, &
       pfft_pencilSize, pfft_procLookup
  implicit none
  include "Flash_mpi.h"
  integer :: ierr, error, i

  do i = 1, pfft_ndim
     deallocate(pfft_procLookup(i) % procInfo, STAT=error)
     if (error /= 0) then
        call Driver_abortFlash("[gr_pfftReconfigFinalise]: " // &
             "Severe error. pfft_procLookup cannot be deallocated!")
     end if
  end do

  call finalise_list(pfft_listFG)  !A list for FLASH grid fragments.
  call finalise_list(pfft_listPG)  !A list for PFFT grid fragments.

  deallocate(pfft_numMsgSendToEachProc)

  call MPI_Type_free(pfft_metaType, ierr)
  call Driver_checkMPIErrorCode(ierr)

  pfft_metaType = MPI_DATATYPE_NULL
  pfft_numFlashNodes = 0; pfft_numPfftNodes = 0
  pfft_srcGridVar = NONEXISTENT; pfft_solnGridVar = NONEXISTENT
  pfft_pencilSize = NONEXISTENT
  nullify(pfft_pfftBuf)

#ifdef DEBUG_PFFT
  call Logfile_close(.true.)
#endif

end subroutine gr_pfftFreeMap
