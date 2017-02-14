!!****f* source/Grid/Grid_renormMassScalars
!!
!! NAME
!!
!!  Grid_renormMassScalars
!!
!!
!! SYNOPSIS
!!
!!  Grid_renormMassScalars(integer(IN)  :: blkLimits(2,MDIM),
!!                         real,pointer :: solnData(:,:,:,:))
!!
!! DESCRIPTION
!!
!!  Renormalize the various Mass Scalar's in groups so they sum to 1.
!!
!! ARGUMENTS
!!
!!  blkLimits - the index limits for internal zones of the block to renormalize
!!  solnData -  Pointer to the block to be renormalized
!!
!!***

subroutine Grid_renormMassScalars(blkLimits,solnData)

  implicit none

#include "constants.h"

  integer, intent(in), dimension(2,MDIM)::blkLimits
  real,pointer :: solnData(:,:,:,:)

  return
end subroutine Grid_renormMassScalars

