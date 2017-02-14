!!****if* source/Grid/GridSolvers/Pfft/MeshReconfiguration/gr_pfftValidateSelectedLevel
!!
!! NAME
!!
!!  gr_pfftValidateSelectedLevel
!!
!! SYNOPSIS
!!  
!!  gr_pfftValidateSelectedLevel(inLevel, outLevel)
!!
!! DESCRIPTION
!!
!!  This subroutine ensures that the selected level can be computed 
!!  in ChooseLevel mode.  In ChooseLevel mode data is restricted / prolonged
!!  from LEAF block level to the chosen level.  As such, we must ensure that 
!!  the LEAF blocks are not resolved too finely relative to the solve level.
!!
!!  Example 1D numerical problem:
!!
!!  Blocks have NXB=8 and lrefine_max=6 and chosen level=2.
!!    At finest resolution the global grid contains 8 * 2**(6-1) = 256 cells.
!!    At level 2 the global grid contains 16 cells.
!!    If a LEAF block exists at level=6 its data must be restricted to level=2.
!!    This does not work because:
!!       Data at level 6: 8 cells.
!!       Data at level 5: Restricted to 4 cells.
!!       Data at level 4: Restricted to 2 cells.
!!       Data at level 3: Restricted to 1 cells.
!!       Data at level 2: Restricted to 0 cells.
!!
!!***

subroutine gr_pfftValidateSelectedLevel(inOutLevel)

use Driver_interface, ONLY : Driver_abortFlash
implicit none
integer, intent(INOUT) :: inOutLevel

call Driver_abortFlash("PFFT Prolong/restrict support not included")

end subroutine gr_pfftValidateSelectedLevel
