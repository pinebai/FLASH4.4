!!****if* source/physics/sourceTerms/Flame/FlameEffects/BurnParametric/fl_effects
!!
!! NAME
!!
!!  fl_effects
!!
!! SYNOPSIS
!!
!!  call fl_effects(real, dimension(:,:,:,:),POINTER_INTENT_IN  :: solndata,
!!                  real,dimension(:,:,:)(in) :: flamdot,
!!                  real(in) :: dt,
!!                  integer(in) :: blockid)
!!
!! DESCRIPTION
!!
!! Dean Townsley 2008
!!
!! save the time derivative representing the actual change in the
!! flam scalar for use by the burn source term
!! TODO should change FLDT to a scratch variable once things have stabilized
!!
!! ARGUMENTS
!!
!!   solndata : 
!!
!!   flamdot : 
!!
!!   dt : 
!!
!!   blockid : ID of block in current processor
!!
!!
!!
!!***



#include "Flash.h"
#include "constants.h"
#include "FortranLangFeatures.fh"
subroutine fl_effects( solnData, flamdot, dt, blockID)

  use Grid_interface, only : Grid_getBlkIndexLimits

  implicit none

  real, dimension(:,:,:,:),POINTER_INTENT_IN  :: solnData
  real,dimension(:,:,:), intent(in)     :: flamdot
  real,intent(in)                       :: dt
  integer, intent(in)                   :: blockID

  integer, dimension(LOW:HIGH,MDIM)     :: blkLimits, blkLimitsGC

  call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)

  ! only need interior cells
  solnData(FLDT_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),   &
                    blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),   &
                    blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) = &
                         flamdot(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),   &
                                 blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),   &
                                 blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))

  return
end subroutine
