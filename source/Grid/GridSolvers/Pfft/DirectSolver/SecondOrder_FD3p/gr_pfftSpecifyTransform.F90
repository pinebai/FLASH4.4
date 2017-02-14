!!****if* source/Grid/GridSolvers/Pfft/DirectSolver/SecondOrder_FD3p/gr_pfftSpecifyTransform
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
!!                 along each axis (e.g. PFFT_REAL2C or PFFT_COMPLEX)
!!
!! NOTES
!!
!! Called by gr_pfftInit.
!!
!!***

subroutine gr_pfftSpecifyTransform(transformType,baseDatType, bcTypes)
  use Grid_interface, ONLY : GRID_PDE_BND_PERIODIC
  use Driver_interface, ONLY : Driver_abortFlash
#include "constants.h"
#include "Flash.h"
#include "Pfft.h"
  implicit none
  integer, dimension(MDIM), intent(OUT) :: transformType
  integer, dimension(0:MDIM), intent(OUT), OPTIONAL :: baseDatType
  integer, dimension(2*MDIM), intent(IN),  OPTIONAL :: bcTypes

  !! In this implementation we are only working with 
  !! periodic boundary conditions
  if (present(bcTypes)) then
     if (bcTypes(1) /= GRID_PDE_BND_PERIODIC .OR. bcTypes(2) /= GRID_PDE_BND_PERIODIC) &
          call Driver_abortFlash("This Poisson solver requires periodic boundaries in the X direction!")
#if NDIM > 1
     if (bcTypes(3) /= GRID_PDE_BND_PERIODIC .OR. bcTypes(4) /= GRID_PDE_BND_PERIODIC) &
          call Driver_abortFlash("This Poisson solver requires periodic boundaries in the Y direction!")
#endif
#if NDIM > 2
     if (bcTypes(5) /= GRID_PDE_BND_PERIODIC .OR. bcTypes(6) /= GRID_PDE_BND_PERIODIC) &
          call Driver_abortFlash("This Poisson solver requires periodic boundaries in the Z direction!")
#endif
  end if

  !! **** THIS IS SPECIFIC TO SecondOrder_FD3p ROUTINE ****
  transformType(IAXIS:KAXIS) = PFFT_COMPLEX

  if (present(baseDatType)) then
     baseDatType(0:MDIM) = PFFT_PCLDATA_COMPLEX
  end if
end subroutine gr_pfftSpecifyTransform
