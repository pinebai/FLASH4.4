!!****f* source/physics/sourceTerms/PrimordialChemistry/PrimordialChemistry_computeDt
!!
!!
!! NAME
!!
!!  PrimordialChemistry_computeDt
!!
!!
!! SYNOPSIS
!!
!!  call PrimordialChemistry_computeDt(integer(IN) :: blockID,
!!                         integer(IN) :: myPE,
!!                         real,pointer :: solnData(:,:,:,:),
!!                         real(INOUT)  :: dt_check,
!!                         integer(INOUT) :: dt_minloc(:) )
!!
!! DESCRIPTION
!!
!!  Computes the timestep limiter.
!!
!!
!! ARGUMENTS
!!
!! blockID      : local block ID
!! solnData     : the physical solution data from grid
!! dt_check     : variable to hold timestep constraint
!! dt_minloc(5) : array to hold limiting zone info: zone indices
!!
!! NOTES
!!
!! stub implementation
!!
!!***

subroutine PrimordialChemistry_computeDt (blockID,  blkLimits, blkLimitsGC, solnData, dt_check, dt_minloc)


#include "constants.h"

        implicit none

        integer, intent(IN) :: blockID
        integer, intent(IN), dimension(2,MDIM) :: blkLimits,blkLimitsGC
        real, pointer :: solnData(:,:,:,:)
        real, INTENT(INOUT) :: dt_check
        integer, INTENT(INOUT) :: dt_minloc(5)

        return
end subroutine PrimordialChemistry_computeDt
