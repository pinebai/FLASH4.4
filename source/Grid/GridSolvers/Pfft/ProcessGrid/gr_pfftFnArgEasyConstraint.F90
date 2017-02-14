!!****if* source/Grid/GridSolvers/Pfft/ProcessGrid/gr_pfftFnArgEasyConstraint
!!
!! NAME
!!
!! gr_pfftFnArgEasyConstraint
!!
!! SYNOPSIS
!!  
!! gr_pfftFnArgEasyConstraint(integer(IN) :: pencilGlobalLen(MDIM), &
!!                            integer(IN) :: totalProcs, &
!!                            integer(IN) :: iProcs, &
!!                            integer(IN) :: jProcs, &
!!                            integer(IN) :: kProcs)
!!
!! DESCRIPTION
!!
!! This is a very easy constraint to satisfy.  We are only interested in 
!! generating a PFFT grid in which we use all available processors.  
!! We should always be able to satisfy this contraint with grids of 
!! [1,1,totalProcs], and [1,totalProcs,1].  Failure to do so is an error! 
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
logical function gr_pfftFnArgEasyConstraint(pencilGlobalLen, totalProcs, &
     iProcs, jProcs, kProcs)
#include "constants.h"
  implicit none
  integer, dimension(MDIM), intent(IN) :: pencilGlobalLen
  integer, intent(IN) :: totalProcs, iProcs, jProcs, kProcs

  gr_pfftFnArgEasyConstraint = .false.
  if (iProcs * jProcs * kProcs == totalProcs) then
     gr_pfftFnArgEasyConstraint = .true.
  end if
end function gr_pfftFnArgEasyConstraint
