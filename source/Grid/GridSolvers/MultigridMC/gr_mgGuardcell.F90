!!****if* source/Grid/GridSolvers/MultigridMC/gr_mgGuardcell
!!
!! NAME
!!
!!  gr_mgGuardcell
!!
!! SYNOPSIS
!!
!!  call gr_mgGuardcell(integer(in) :: mype2,
!!                      integer(in) :: ivar,
!!                      integer(in) :: nlayers,
!!                      real (in) :: simtime,
!!                      integer(in) :: idiag,
!!                      integer(in) :: idir)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   mype2 : 
!!
!!   ivar : 
!!
!!   nlayers : 
!!
!!   simtime : 
!!
!!   idiag : 
!!
!!   idir : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

!*******************************************************************************

! Routine:      mg_guardcell

! Description:  Version of amr_guardcell() for the multigrid solver.
!               The only differences are (1) it only operates on the
!               work array and (2) it calls the special routine
!               mg_set_ext_bndry() to handle exterior guard cells.

!               See documentation for amr_guardcell() for more information.


subroutine gr_mgGuardcell(mype2, ivar, nlayers, simtime, idiag, idir)

!===============================================================================

use gr_mgData, only :  gr_mgDiffOpDiscretize

use Driver_interface, ONLY : Driver_abortFlash
use paramesh_interfaces, only: amr_guardcell
use paramesh_dimensions, only : nguard_work

use physicaldata, only : diagonals

use workspace, only : interp_mask_work,interp_mask_work_res

use Grid_interface, ONLY : Grid_fillGuardCells

implicit none

!include 'mpif.h'
#include "Flash.h"
#include "constants.h"

integer, intent(in) :: mype2, nlayers, idiag, idir
integer, intent(in) :: ivar
real   , intent(in) :: simtime
!===============================================================================

if (nlayers > nguard_work) &
  call Driver_abortFlash("mg_guardcell: you must set nlayers <= nguard_work")


select case(gr_mgDiffOpDiscretize)
case(2)
!interp_mask_work(:) = 2
!interp_mask_work_res(:) = 2
diagonals = .false.
call amr_guardcell(mype2,2,nlayers,nlayersx=1,nlayersy=1,nlayersz=1,maxNodetype_gcWanted=1)
diagonals = .true.
!call Grid_fillGuardCells( WORK, ALLDIR, minLayers=1)
!interp_mask_work(:) = 2
!interp_mask_work_res(:) = 2
case(4)
!interp_mask_work(:) = 4
!interp_mask_work_res(:) = 2
diagonals = .false.
call amr_guardcell(mype2,2,nlayers,nlayersx=2,nlayersy=2,nlayersz=2)
diagonals = .true.
!call Grid_fillGuardCells( WORK, ALLDIR, minLayers=3)
!interp_mask_work(:) = 2
!interp_mask_work_res(:) = 2
end select

!===============================================================================

return
end
