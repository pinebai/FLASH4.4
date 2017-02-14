!!****if* source/Grid/GridSolvers/Pfft/ProcessGrid/gr_pfftFnArgMediumConstraint
!!
!! NAME
!!
!! gr_pfftFnArgMediumConstraint
!!
!! SYNOPSIS
!!  
!! gr_pfftFnArgMediumConstraint(integer(IN) :: pencilGlobalLen(MDIM), &
!!                              integer(IN) :: totalProcs, &
!!                              integer(IN) :: iProcs, &
!!                              integer(IN) :: jProcs, &
!!                              integer(IN) :: kProcs)
!!
!! DESCRIPTION
!!
!! This is a difficult contraint to satisfy because the number of global 
!! grid points in a single dimension must be divisible by the number of 
!! PFFT processors in the same dimension.
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
logical function gr_pfftFnArgMediumConstraint(pencilGlobalLen, totalProcs, &
     iProcs, jProcs, kProcs)
#include "constants.h"
  implicit none
  integer, dimension(MDIM), intent(IN) :: pencilGlobalLen
  integer, intent(IN) :: totalProcs, iProcs, jProcs, kProcs

  gr_pfftFnArgMediumConstraint = .false.
  if (iProcs * jProcs * kProcs == totalProcs) then
     if ( (mod(pencilGlobalLen(IAXIS), iProcs) == 0).and.&
          (mod(pencilGlobalLen(JAXIS), jProcs) == 0).and.&
          (mod(pencilGlobalLen(KAXIS), kProcs) == 0) ) then
        gr_pfftFnArgMediumConstraint = .true.
     end if
  end if
end function gr_pfftFnArgMediumConstraint
