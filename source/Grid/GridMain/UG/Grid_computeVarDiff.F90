!!****if* source/Grid/GridMain/UG/Grid_computeVarDiff
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

!!REORDER(5): unk

subroutine Grid_computeVarDiff(level, gr_iRefSoln, gr_iSoln, ires)

!==============================================================================
#include "Flash.h"
#include "constants.h"

  use Grid_interface, ONLY : Grid_getBlkIndexLimits
  use physicaldata, ONLY : unk

  use Timers_interface, ONLY : Timers_start, Timers_stop

  implicit none

  integer, intent(in)          :: level, gr_iRefSoln, gr_iSoln, ires
  
  integer, dimension(2,MDIM)   :: blkLimitsGC, blkLimits
  integer                      :: i, j, k

  
  !=======================================================================

  call Timers_start("Grid_computeVarDiff")
  
  
  
  call Grid_getBlkIndexLimits(1,blkLimits,blkLimitsGC)

  if (level==-1 .OR. level == 1) then

     do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)  ! working on interior only
        do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
              unk(ires,i,j,k,1) = unk(gr_iSoln, i,j,k,1) - unk(gr_iRefSoln, i,j,k,1)
           enddo
        enddo
     enddo
!!  else
     ! DEV: Issue warning?
  endif


  call Timers_stop("Grid_computeVarDiff")

  !=======================================================================
  
  return
end subroutine Grid_computeVarDiff
