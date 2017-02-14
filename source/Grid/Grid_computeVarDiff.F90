!!****f* source/Grid/Grid_computeVarDiff
!!
!! NAME
!!  Grid_computeVarDiff
!!
!! SYNOPSIS
!!  
!!  call Grid_computeVarDiff(integer(in) :: level,
!!                integer(in) :: gr_iRefSoln,
!!                integer(in) :: gr_iSoln,
!!                integer(in) :: ires)
!!
!! DESCRIPTION
!!  
!!  Compute the current difference of two variables.
!!  If level is -1, only the leaf block residuals are calculated;
!!  otherwise, residuals are computed for all blocks of the given level.
!!
!!
!! ARGUMENTS
!!
!!  level        - the level of blocks to take the residual for, or -1
!!                 for all LEAF blocks.
!!  gr_iRefSoln  - the "reference solution" variable
!!  gr_iSoln     - "the solution" variable
!!  ires         - the residual variable
!!
!! RESULT
!!
!!  The difference between the reference solution and solution in the blocks
!!  (either selected by level, or all LEAF blocks)
!!  is placed in the variable ires.
!!
!!***

subroutine Grid_computeVarDiff(level, gr_iRefSoln, gr_iSoln, ires)

  implicit none

  integer, intent(in)          :: level, gr_iRefSoln, gr_iSoln, ires

end subroutine Grid_computeVarDiff
