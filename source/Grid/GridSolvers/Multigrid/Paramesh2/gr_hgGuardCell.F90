!!****if* source/Grid/GridSolvers/Multigrid/Paramesh2/gr_hgGuardCell
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
!!  paramesh4 Grid_fillGuardcells.
!!
!!***

Subroutine gr_hgGuardCell(myPE, nlayers, extrap)

  use gr_hgData, ONLY: gr_hgNowActive, gr_hgCurrentGcReq_extrap
  use workspace, ONLY: nguard_work
  use Grid_data, ONLY : gr_intpolStencilWidth

  use paramesh_interfaces, ONLY: amr_restrict
  use Timers_interface, ONLY: Timers_start, Timers_stop
  use Driver_interface, ONLY: Driver_abortFlash
  use gr_hgInterface, ONLY: gr_hgSetExtBoundary
  
  implicit none

#include "constants.h"
#include "Flash.h"

  integer, INTENT(in) :: myPE, nlayers
  logical, intent(in) :: extrap
  integer :: nlayers_transverse_parents, nlayersAllLevels, nlayersForInterp
  integer :: idiag, idir

!=====================================================================
  idiag = 1                     ! hardwired - fill diagonals.
  idir = 0                      ! hardwired - fill in all directions.

  call Timers_start("gr_hgGuardCell")
  gr_hgNowActive = .TRUE.
  gr_hgCurrentGcReq_extrap = extrap
  
if (nlayers > NGUARD_WORK) &
  call Driver_abortFlash("gr_hgGuardCell: you must set nlayers <= NGUARD_work")

! Put recognizably bad data in corners to permit recognition of
! corner guard cells diagonally opposite coarser blocks.

if (idiag == 1) call amr_mark_edges(mype, 2)

! Restrict data from leaf blocks to their parents. This will enable the
! parent blocks to provide guard cell data to their neighbors in cases
! where the leaf blocks have coarser neighbors.

call amr_restrict(mype, 2, 0)

!!  nlayers must be an even number for monotonic interpolation,
!!  and may be changed to ensure that parent blocks have enough guard cells for interpolation.
!!  The fixes here are similar to the logic in Grid_setGcFillNLayer, which is
!!  used with paramesh4.

nlayers_transverse_parents = ((nlayers + 1) / 2) + gr_intpolStencilWidth
nlayersAllLevels = max(nlayers,nlayers_transverse_parents)
  
  ! Blocks provide guard cell data to all their neighbors which share their
  ! level of refinement.
  
  !call timer_start("guardcell srl")
  call amr_guardcell_srl(mype, 2, nlayersAllLevels, idiag, idir)
  !call timer_stop("guardcell srl")

  ! Apply boundary conditions by filling guard cells at external boundaries.
  ! This first call to gr_hgSetExtBoundary is done only to prepare the
  ! input for amr_guardcell_c_to_f interpolation.
  call gr_hgSetExtBoundary(idiag, idir, extrap=.TRUE.)

  ! Guard cell data is sent to any neighbor blocks with finer resolution.
  nlayersForInterp = nlayers
#ifdef GRID_GC_LAYERS_ALWAYS_EVEN
     nlayersForInterp = ((nlayersForInterp + 1) / 2)  * 2
#endif
#ifdef DEBUG_SOLVERS
   print*,'NL,NLaL,NLfI:',nlayers,nlayersAllLevels,nlayersForInterp
#endif
  !call timer_start("guardcell c_to_f")
  call amr_guardcell_c_to_f(mype, 2, nlayersForInterp, idiag, idir)
  !call timer_stop("guardcell c_to_f")


  ! Apply boundary conditions by filling guard cells at external boundaries.
  ! This second call to gr_hgSetExtBoundary is done to fill the guard cells
  ! in the way the caller actually requests.
  ! (This call may be unnecessary if extrap is .TRUE., since in that case
  ! guard cells have already been filled by the first call above! - KW)
  call gr_hgSetExtBoundary(1, 0, extrap)

  call Timers_stop("gr_hgGuardCell")
  gr_hgNowActive = .FALSE.

  return
End subroutine gr_hgGuardCell

