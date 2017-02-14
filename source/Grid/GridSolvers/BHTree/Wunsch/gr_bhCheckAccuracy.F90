!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhCheckAccuracy
!!
!! NAME
!!
!!  gr_bhCheckAccuracy
!!
!!
!! SYNOPSIS
!!
!!  logical res = gr_bhCheckAccuracy(
!!                          integer(in)    :: blockno,
!!                          integer(in)    :: point(MDIM),
!!                          integer(in)    :: blkLimits(2,MDIM),
!!                          real,pointer   :: solnData(:,:,:,:)
!!        )
!!
!! DESCRIPTION
!!
!!  Used by gr_bhTreeWalkBlockUnified to evaluate whether quantities within 
!!  a given block calculated by tree solver are accurate enough and do not 
!!  need to be recalculated.
!!
!! ARGUMENTS
!!
!!  blockno     : number of block into which the target cell belongs
!!  point       : indeces of the target cell in the block
!!  blkLimits   : limits of indeces in the block
!!  solnData    : solution data from the grid
!!
!! RESULT
!!
!!  Returns TRUE if quantities within a given block calculated by tree solver 
!!  are accurate enough and do not need to be recalculated.
!!
!!***

logical function gr_bhCheckAccuracy(blockno, point, blkLimits, solnData)
  use Gravity_interface, ONLY : Gravity_bhCheckAccuracy
  use TreeRay_interface, ONLY : TreeRay_bhCheckAccuracy
  use gr_bhData, ONLY : gr_bhUseGravity, gr_bhUseTreeRay
  implicit none
#include "FortranLangFeatures.fh"
#include "constants.h"
  integer, intent(IN) :: blockno
  integer, dimension(MDIM), intent(IN) :: point
  integer, dimension(2,MDIM), intent(IN)   :: blkLimits
  real, DIMENSION(:,:,:,:), POINTER_INTENT_IN :: solnData
  logical :: res

  res = .true.

  ! CheckAccuracy of physical modules
  ! res is tested because fortran does not ensure short-circuit .AND. evaluation
  if (gr_bhUseGravity .and. res) &
  & res = res .and. Gravity_bhCheckAccuracy(blockno, point, blkLimits, solnData)
  if (gr_bhUseTreeRay .and. res) &
  & res = res .and. TreeRay_bhCheckAccuracy(blockno, point, blkLimits, solnData)

  gr_bhCheckAccuracy = res

end function gr_bhCheckAccuracy

