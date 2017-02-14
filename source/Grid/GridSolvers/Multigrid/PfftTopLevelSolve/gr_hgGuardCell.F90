!!****if* source/Grid/GridSolvers/Multigrid/PfftTopLevelSolve/gr_hgGuardCell
!!
!! NAME
!!  gr_hgGuardCell
!!
!! SYNOPSIS
!!
!!  gr_hgGuardCell(integer, intent(in) :: myPE,
!!                 integer, intent(in) :: nlayers,
!!                 logical, intent(in) :: extrap)
!!
!! DESCRIPTION
!!
!!  This routine fills the guardcells in the work variable, and 
!!  sets the external boundaries appropriately.
!!
!! ARGUMENTS
!!
!!  myPE    - the current processor
!!  nlayers - the number of layers to fill
!!  extrap  - extrapolate the cell values into the outer boundaries
!!              .false. => fill the boundary as if it goes to zero
!!              .true.  => extrapolate the interior cell values outwards
!!
!! NOTES
!!
!!  calling amr_guardcell with inadequate numbers of layers to fill creates 
!!  problems with monotonic interpolation.  The fix below echoes that of
!!  Paramesh3/Grid_fillGuardcells
!!
!!***

Subroutine gr_hgGuardCell(myPE, nlayers, extrap)

  use workspace, ONLY: interp_mask_work
  use gr_hgData, ONLY: gr_hgNowActive, gr_hgCurrentGcReq_extrap

  use paramesh_interfaces, ONLY: amr_guardcell
  use Timers_interface, ONLY: Timers_start, Timers_stop
  use gr_hgInterface, ONLY: gr_hgSetExtBoundary
  use gr_interface, ONLY : gr_setGcFillNLayers
  
  implicit none

#include "constants.h"
#include "Flash.h"

  integer, INTENT(in) :: myPE, nlayers
  logical, intent(in) :: extrap
  integer :: nlayersx, nlayersy, nlayersz

!=====================================================================

  call Timers_start("gr_hgGuardCell")
  gr_hgNowActive = .TRUE.
  gr_hgCurrentGcReq_extrap = extrap
  
!  Klaus says that nlayers must be an even number in amr_guardCell for monotonic interpolation,
!  and may be changed to ensure that parent blocks have enough guard cells for interpolation.
!  The fix here follows that of Grid_fillGuardCells
!  call gr_setGcFillNLayers(nlayersx,nlayersy,nlayersz, ALLDIR, NGUARD, nlayers)
  nlayersx = NGUARD 
  nlayersy = NGUARD
  nlayersz = NGUARD

  call amr_guardcell(myPE, 2, NGUARD, nlayersx, nlayersy, nlayersz)

! Should not be needed - should verify that results do not depend on it - KW
  call gr_hgSetExtBoundary(1, 0, extrap)

#ifdef DEBUG_GRID
  print*,'After amr_guardcell, call gr_hgSetExtBoundary(1, 0,', extrap,')'  ! Within DEBUG
#endif 

  call Timers_stop("gr_hgGuardCell")
  gr_hgNowActive = .FALSE.

  return
End subroutine gr_hgGuardCell

