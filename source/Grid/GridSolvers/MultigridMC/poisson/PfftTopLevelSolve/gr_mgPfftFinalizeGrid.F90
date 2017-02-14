!!****if* source/Grid/GridSolvers/MultigridMC/poisson/PfftTopLevelSolve/gr_mgPfftFinalizeGrid
!!
!! NAME
!!
!!  gr_mgPfftFinalizeGrid
!!
!! SYNOPSIS
!!
!!  gr_mgPfftFinalizeGrid()
!!
!! DESCRIPTION
!!
!!  This routine cleans up the data necessary for the Multigrid PFFT extensions.
!!
!!***

subroutine gr_mgPfftFinalizeGrid()

  use gr_mgPfftData, ONLY : gr_mgPfftInArray, gr_mgPfftOutArray, gr_mgPfftTranArray, &
       gr_mgPfftLastMappedLevel
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_pfftFinalize

#include "constants.h"
  implicit none
  integer :: error

  if (gr_mgPfftLastMappedLevel /= NONEXISTENT) then
     call Grid_pfftFinalize()

     !Only allocated if processor is part of PFFT communicator.
     if (allocated(gr_mgPfftInArray)) then
        deallocate(gr_mgPfftInArray, STAT=error)
        if (error /= 0) then
           call Driver_abortFlash &
                ("[gr_mgPfftFinalizeGrid]: gr_mgPfftInArray cannot be deallocated!")
        end if
     end if

     if (allocated(gr_mgPfftOutArray)) then
        deallocate(gr_mgPfftOutArray, STAT=error)
        if (error /= 0) then
           call Driver_abortFlash &
                ("[gr_mgPfftFinalizeGrid]: gr_mgPfftOutArray cannot be deallocated!")
        end if
     end if
  end if

end subroutine gr_mgPfftFinalizeGrid
