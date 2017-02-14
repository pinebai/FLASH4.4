!!****if* source/Grid/GridSolvers/Pfft/ImplicitDiffusionSolver/gr_pfftSpecifyTransform
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

subroutine gr_pfftSpecifyTransform(transformType,baseDatType,bcTypes)
#include "constants.h"
#include "Pfft.h"
  implicit none
  integer, dimension(MDIM), intent(OUT) :: transformType
  integer, dimension(0:MDIM), intent(OUT), OPTIONAL :: baseDatType
  integer, dimension(2*MDIM), intent(IN),  OPTIONAL :: bcTypes

  !! This is some trickery to convince the Pfft machinery that
  !! we do not have complex data, neither on input nor for output
  !! nor at any point in between.
  transformType(IAXIS) = PFFT_SIN
  transformType(JAXIS:KAXIS) = PFFT_SIN

  if (present(baseDatType)) then
     baseDatType(:) = PFFT_PCLDATA_REAL
  end if
end subroutine gr_pfftSpecifyTransform
