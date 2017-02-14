!!****if* source/Grid/GridSolvers/Pfft/ProcessGrid/gr_pfftFnArgHardConstraint
!!
!! NAME
!!
!! gr_pfftFnArgHardConstraint
!!
!! SYNOPSIS
!!  
!! gr_pfftFnArgHardConstraint(integer(IN) :: pencilGlobalLen(MDIM), &
!!                            integer(IN) :: totalProcs, &
!!                            integer(IN) :: iProcs, &
!!                            integer(IN) :: jProcs, &
!!                            integer(IN) :: kProcs)
!!
!! DESCRIPTION
!!
!! This is the hardest set of contraints to satisfy: The global grid points 
!! must be divisible by the number of PFFT processors in the same dimension.  
!! We also require that the global grid points must be divisible by the number 
!! of FLASH grid processors in the same dimension.
!!
!! ARGUMENTS
!!
!! pencilGlobalLen - Array containing no. cells in global pencil grid.
!! totalProcs - The total number of processors we must use.
!! iProcs - A candidate value for the number of IAXIS processors.
!! jProcs - A candidate value for the number of JAXIS processors.
!! kProcs - A candidate value for the number of KAXIS processors.
!!
!!***
logical function gr_pfftFnArgHardConstraint(pencilGlobalLen, totalProcs, &
     iProcs, jProcs, kProcs)
#include "constants.h"
#include "Flash.h"

#if defined(FLASH_GRID_UG)
  use Grid_data, ONLY : gr_axisNumProcs
#elif defined(FLASH_GRID_PARAMESH)
  use Grid_data, ONLY : gr_ilo, gr_ihi, gr_jlo, gr_jhi, gr_klo, gr_khi
#endif
  use Driver_interface, ONLY : Driver_abortFlash
  implicit none
  integer, dimension(MDIM), intent(IN) :: pencilGlobalLen
  integer, intent(IN) :: totalProcs, iProcs, jProcs, kProcs
  integer :: flashIBlocks, flashJBlocks, flashKBlocks

#if defined(FLASH_GRID_UG)
  flashIBlocks = gr_axisNumProcs(IAXIS)
  flashJBlocks = gr_axisNumProcs(JAXIS)
  flashKBlocks = gr_axisNumProcs(KAXIS)
#elif defined(FLASH_GRID_PARAMESH)
  !This is the number of blocks across the computational domain
  !at the chosen refinement level.  Remember we adjusted 
  !pfft_globalLen earlier to be consistent with gr_oneRefLev.
  flashIBlocks = pencilGlobalLen(IAXIS) / (gr_ihi - gr_ilo + 1)
  flashJBlocks = pencilGlobalLen(JAXIS) / (gr_jhi - gr_jlo + 1)
  flashKBlocks = pencilGlobalLen(KAXIS) / (gr_khi - gr_klo + 1)
#else
  call Driver_abortFlash("[gr_pfftFnArgHardConstraint]: "//&
       "Situation not yet considered")
#endif

  gr_pfftFnArgHardConstraint = .false.
  if (iProcs * jProcs * kProcs == totalProcs) then
     if ( (mod(pencilGlobalLen(IAXIS), iProcs) == 0).and.&
          (mod(pencilGlobalLen(JAXIS), jProcs) == 0).and.&
          (mod(pencilGlobalLen(KAXIS), kProcs) == 0).and.&
          (mod(pencilGlobalLen(IAXIS), flashIBlocks) == 0).and.&
          (mod(pencilGlobalLen(JAXIS), flashJBlocks) == 0).and.&
          (mod(pencilGlobalLen(KAXIS), flashKBlocks) == 0) ) then
        gr_pfftFnArgHardConstraint = .true.
     end if
  end if
end function gr_pfftFnArgHardConstraint
