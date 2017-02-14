!!****if* source/physics/sourceTerms/Flame/FlameSpeed/Constant/TFI/fl_flameSpeed
!!
!! NAME
!!
!!  fl_flameSpeed
!!
!! SYNOPSIS
!!
!!  call fl_flameSpeed(real, dimension(:,:,:,:),POINTER_INTENT_IN  :: solndata,
!!                     real, dimension(:,:,:)(out) :: flamespeed,
!!                     integer(in) :: blockid,
!!                     integer(in) :: nlayers)
!!
!! DESCRIPTION
!!
!! Dean Townsley 2008
!!
!! ARGUMENTS
!!
!!   solndata : solution data
!!
!!   flamespeed : flame speed
!!
!!   blockid : ID of block in current processor
!!
!!   nlayers : number of layers
!!
!!
!!
!!***


subroutine fl_flameSpeed( solnData, flamespeed, blockID, nlayers)

#include "Flash.h"
#include "constants.h"
#include "FortranLangFeatures.fh"

  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getDeltas
  use fl_fsData, only : fl_fsConstFlameSpeed, fl_fsConstFlameWidth, &
        fl_fsUseTFI
  use fl_fsTFIInterface, ONLY : fl_fsTFIFlameSpeedBlock
  implicit none
  real, dimension(:,:,:,:),POINTER_INTENT_IN :: solnData
  real, dimension(:,:,:),intent(out) :: flamespeed
  integer, intent(in) :: blockID, nlayers

  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC, compLimits
  real, dimension(MDIM)             :: deltas

#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_ILO_GC:GRID_IHI_GC,&
                  GRID_JLO_GC:GRID_JHI_GC,&
                  GRID_KLO_GC:GRID_KHI_GC) :: flamewidth
#else
  real, allocatable, dimension(:,:,:) :: flamewidth
  integer :: sizeI, sizeJ, sizeK, istat
#endif

  real :: dx

  if (.not. fl_fsUseTFI) then
     flamespeed(:,:,:) = fl_fsConstFlameSpeed
#ifdef FSPD_VAR
     solndata(FSPD_VAR,:,:,:) = fl_fsConstFlameSpeed
#endif
     return
  endif

  call Grid_getDeltas(blockID, deltas)
  ! assume square grid
  dx = deltas(IAXIS)

  call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)

 ! make an array of limits over which computation is being performed
  compLimits(LOW,IAXIS) = blkLimits(LOW,IAXIS)-nlayers
  compLimits(HIGH,IAXIS) = blkLimits(HIGH,IAXIS)+nlayers
  compLimits(LOW,JAXIS) = blkLimits(LOW,JAXIS)-K2D*nlayers
  compLimits(HIGH,JAXIS) = blkLimits(HIGH,JAXIS)+K2D*nlayers
  compLimits(LOW,KAXIS) = blkLimits(LOW,KAXIS)-K3D*nlayers
  compLimits(HIGH,KAXIS) = blkLimits(HIGH,KAXIS)+K3D*nlayers

#ifndef FIXEDBLOCKSIZE
  sizeI=blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  sizeJ=blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
  sizeK=blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1
  allocate(flamewidth(sizeI,sizeJ,sizeK),STAT=istat)
  if (istat /= 0) call Driver_abortFlash("Cannot allocate flamewidth in fl_flameSpeed")
#endif

  flamespeed( compLimits(LOW,IAXIS):compLimits(HIGH,IAXIS), &
              compLimits(LOW,JAXIS):compLimits(HIGH,JAXIS), &
              compLimits(LOW,KAXIS):compLimits(HIGH,KAXIS)) &
              = fl_fsConstFlameSpeed
  flamewidth( compLimits(LOW,IAXIS):compLimits(HIGH,IAXIS), &
              compLimits(LOW,JAXIS):compLimits(HIGH,JAXIS), &
              compLimits(LOW,KAXIS):compLimits(HIGH,KAXIS)) &
              = fl_fsConstFlameWidth

  call fl_fsTFIFlameSpeedBlock(solnData, flamespeed, flamewidth, dx, &
          compLimits)

#ifdef FSPD_VAR
     solndata(FSPD_VAR, compLimits(LOW,IAXIS):compLimits(HIGH,IAXIS), &
                        compLimits(LOW,JAXIS):compLimits(HIGH,JAXIS), &
                        compLimits(LOW,KAXIS):compLimits(HIGH,KAXIS)) &
                  = flamespeed(compLimits(LOW,IAXIS):compLimits(HIGH,IAXIS), &
                               compLimits(LOW,JAXIS):compLimits(HIGH,JAXIS), &
                               compLimits(LOW,KAXIS):compLimits(HIGH,KAXIS))
#endif

#ifndef FIXEDBLOCKSIZE
  deallocate(flamewidth)
#endif

  return
end subroutine
