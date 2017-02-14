!!****f* source/physics/Gravity/Gravity_bhCheckAccuracy
!!
!! NAME
!!
!!  Gravity_bhCheckAccuracy
!!
!!
!! SYNOPSIS
!!
!!  logical res = Gravity_bhCheckAccuracy(
!!                          integer(in)    :: blockno,
!!                          integer(in)    :: point(MDIM),
!!                          integer(in)    :: blkLimits(2,MDIM),
!!                          real,pointer   :: solnData(:,:,:,:)
!!        )
!!
!! DESCRIPTION
!!
!!  Evaluates whether gravitational acceleration/potential within 
!!  a given block calculated by tree solver is accurate enough and 
!!  does not need to be recalculated.
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
!!  Returns TRUE if the gravitational potential/aceleration within 
!!  a given block calculated by tree solver is accurate enough and 
!!  does not need to be recalculated.
!!
!!***

logical function Gravity_bhCheckAccuracy(blockno, point, blkLimits, solnData)
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
  integer, intent(IN) :: blockno
  integer, dimension(MDIM), intent(IN) :: point
  integer, dimension(2,MDIM), intent(IN)   :: blkLimits
  real, DIMENSION(:,:,:,:), POINTER :: solnData

  Gravity_bhCheckAccuracy = .true.
  return
end function Gravity_bhCheckAccuracy


