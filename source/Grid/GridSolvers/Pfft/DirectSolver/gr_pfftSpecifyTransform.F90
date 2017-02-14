!!****if* source/Grid/GridSolvers/Pfft/DirectSolver/gr_pfftSpecifyTransform
!!
!! NAME
!!
!! gr_pfftSpecifyTransform
!!
!! SYNOPSIS
!!  
!! gr_pfftSpecifyTransform(integer(OUT) :: transformType(MDIM))
!! 
!! DESCRIPTION
!!
!! Initialises the transform type for the particular solver.
!!
!! ARGUMENTS
!!  
!! transformType - Array containing the type of transform required 
!!                 along each axis (e.g. PFFT_REAL or PFFT_COMPLEX)
!!
!! NOTES
!!
!! Called by gr_pfftInit.
!!
!!***

subroutine gr_pfftSpecifyTransform(transformType,baseDatType, bcTypes)
  use Grid_interface, ONLY : GRID_PDE_BND_PERIODIC, GRID_PDE_BND_NEUMANN
  use Driver_interface, ONLY : Driver_abortFlash
#include "constants.h"
#include "Flash.h"
#include "Pfft.h"
  implicit none
  integer, dimension(MDIM), intent(OUT) :: transformType
  integer, dimension(0:MDIM), intent(OUT), OPTIONAL :: baseDatType
  integer, dimension(2*MDIM), intent(IN),  OPTIONAL :: bcTypes

  !! **** THIS IS SPECIFIC TO DIRECTSOLVER ROUTINE ****

  if (present(bcTypes)) then
     if (bcTypes(1) /= GRID_PDE_BND_PERIODIC .OR. bcTypes(2) /= GRID_PDE_BND_PERIODIC) &
          call Driver_abortFlash("This Poisson solver requires periodic boundaries in the X direction!")
#if NDIM > 1
     if (bcTypes(3) /= bcTypes(4)) &
          call Driver_abortFlash("This Poisson solver requires the same type of boundaries up and down for the Y direction!")
     if (bcTypes(3) /= GRID_PDE_BND_PERIODIC .AND. bcTypes(3) /= GRID_PDE_BND_NEUMANN) &
          call Driver_abortFlash("This Poisson solver requires periodic or Neumann boundaries in the Y direction!")
#endif
#if NDIM > 2
     if (bcTypes(5) /= bcTypes(6)) &
          call Driver_abortFlash("This Poisson solver requires the same type of boundaries left and right for the Z direction!")
     if (bcTypes(5) /= GRID_PDE_BND_PERIODIC .AND. bcTypes(5) /= GRID_PDE_BND_NEUMANN) &
          call Driver_abortFlash("This Poisson solver requires periodic or Neumann boundaries in the Z direction!")
#endif
  end if

  transformType(IAXIS:KAXIS) = PFFT_REAL

  if (present(baseDatType)) then
     baseDatType(0:MDIM) = PFFT_PCLDATA_REAL
  end if
end subroutine gr_pfftSpecifyTransform
