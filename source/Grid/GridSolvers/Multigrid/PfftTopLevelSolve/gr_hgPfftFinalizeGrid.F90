!!****if* source/Grid/GridSolvers/Multigrid/PfftTopLevelSolve/gr_hgPfftFinalizeGrid
!!
!! NAME
!!
!!  gr_hgPfftFinalizeGrid
!!
!! SYNOPSIS
!!
!!  gr_hgPfftFinalizeGrid()
!!
!! DESCRIPTION
!!
!!  This routine cleans up the data necessary for the Multigrid PFFT extensions.
!!
!!***

subroutine gr_hgPfftFinalizeGrid()

  use gr_hgPfftData, ONLY : gr_hgPfftInArray, gr_hgPfftOutArray, gr_hgPfftTranArray, &
       gr_hgPfftLastMappedLevel
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_pfftFinalize

#include "constants.h"
  implicit none
  integer :: error

  if (gr_hgPfftLastMappedLevel /= NONEXISTENT) then
     call Grid_pfftFinalize()

     !Only allocated if processor is part of PFFT communicator.
     if (allocated(gr_hgPfftInArray)) then
        deallocate(gr_hgPfftInArray, STAT=error)
        if (error /= 0) then
           call Driver_abortFlash &
                ("[gr_hgPfftFinalizeGrid]: gr_hgPfftInArray cannot be deallocated!")
        end if
     end if

     if (allocated(gr_hgPfftOutArray)) then
        deallocate(gr_hgPfftOutArray, STAT=error)
        if (error /= 0) then
           call Driver_abortFlash &
                ("[gr_hgPfftFinalizeGrid]: gr_hgPfftOutArray cannot be deallocated!")
        end if
     end if
  end if

end subroutine gr_hgPfftFinalizeGrid
