!!****if* source/physics/sourceTerms/Flame/FlameSpeed/LaminarOnly/fl_flameSpeed
!!
!! NAME
!!
!!  fl_flameSpeed
!!
!! SYNOPSIS
!!
!!  call fl_flameSpeed(real,dimension(:,:,:,:),pointer  :: solndata,
!!                     real,dimension(:,:,:)(out) :: flamespeed,
!!                     integer(in) :: blockid,
!!                     integer(in) :: nlayers)
!!
!! DESCRIPTION
!!  Aaron Jackson, Dean Townsley, Alan Calder 2008
!!  current implementation ignores nlayers
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
!!***


subroutine fl_flameSpeed(solnData, flamespeed, blockID, nlayers)

#include "constants.h"
#include "Flash.h"
  
  use fl_fsData, ONLY : fl_fsUseConstFlameSpeed, fl_fsConstFlameSpeed, &
         fl_fsConstFlameWidth, fl_fsUseTFI
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getDeltas
  use fl_fsLaminarInterface, ONLY : fl_fsLaminarFlameSpeedBlock
  use fl_fsTFIInterface, ONLY : fl_fsTFIFlameSpeedBlock
  implicit none
  
  real,dimension(:,:,:,:),pointer :: solnData
  real,dimension(:,:,:),intent(out) :: flamespeed
  integer, intent(in) :: blockID, nlayers
  
  integer, dimension(LOW:HIGH,MDIM):: blkLimits, blkLimitsGC, compLimits
  real, dimension(MDIM)            :: deltas
#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_ILO_GC:GRID_IHI_GC,&
                  GRID_JLO_GC:GRID_JHI_GC,&
                  GRID_KLO_GC:GRID_KHI_GC) :: dens, flamewidth
#else
  real, allocatable, dimension(:,:,:) :: dens, flamewidth
  integer :: sizeX, sizeY, sizeZ, istat
#endif

  real :: dx

!.. jump out for constant flame speed
  if (fl_fsUseConstFlameSpeed .and. (.not. fl_fsUseTFI)) then
     flamespeed(:,:,:) = fl_fsConstFlameSpeed
#ifdef FSPD_VAR
     solndata(FSPD_VAR,:,:,:) = fl_fsConstFlameSpeed
#endif
     return
  endif

  call Grid_getDeltas(blockID, deltas)
  ! assume square grid
  dx = deltas(IAXIS)

  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

#ifndef FIXEDBLOCKSIZE
  sizeX=blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  sizeY=blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
  sizeZ=blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1
  allocate(flamewidth(sizeX,sizeY,sizeZ),STAT=istat)
  if (istat /= 0) call Driver_abortFlash("Cannot allocate flamewidth in fl_flameSpeed")
#endif

  ! make an array of limits over which computation is being performed
  compLimits(LOW,IAXIS) = blkLimits(LOW,IAXIS)-nlayers
  compLimits(HIGH,IAXIS) = blkLimits(HIGH,IAXIS)+nlayers
  compLimits(LOW,JAXIS) = blkLimits(LOW,JAXIS)-K2D*nlayers
  compLimits(HIGH,JAXIS) = blkLimits(HIGH,JAXIS)+K2D*nlayers
  compLimits(LOW,KAXIS) = blkLimits(LOW,KAXIS)-K3D*nlayers
  compLimits(HIGH,KAXIS) = blkLimits(HIGH,KAXIS)+K3D*nlayers


!.. jump out for constant flame speed with TFI
  if (fl_fsUseConstFlameSpeed) then
     flamespeed(:,:,:) = fl_fsConstFlameSpeed
     flamewidth(:,:,:) = fl_fsConstFlameWidth
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
  endif

#ifndef FIXEDBLOCKSIZE
  allocate(dens(sizeX,sizeY,sizeZ),STAT=istat)
  if (istat /= 0) call Driver_abortFlash("Cannot allocate dens in fl_flameSpeed")
#endif

  ! this returns an estimate of unburned density at this pressure
  call fl_fsLaminarFlameSpeedBlock(solnData, flamespeed, flamewidth, dens,  &
                                   compLimits)
  if (fl_fsUseTFI) &
    call fl_fsTFIFlameSpeedBlock(solnData, flamespeed, flamewidth, dx, &
      compLimits)

#ifdef FSPD_VAR
  !.. save flame speed for debugging
  solndata(FSPD_VAR, compLimits(LOW,IAXIS):compLimits(HIGH,IAXIS), &
                     compLimits(LOW,JAXIS):compLimits(HIGH,JAXIS), &
                     compLimits(LOW,KAXIS):compLimits(HIGH,KAXIS)) &
                 = flamespeed(compLimits(LOW,IAXIS):compLimits(HIGH,IAXIS), &
                              compLimits(LOW,JAXIS):compLimits(HIGH,JAXIS), &
                              compLimits(LOW,KAXIS):compLimits(HIGH,KAXIS))
#endif

#ifndef FIXEDBLOCKSIZE
  deallocate(dens)
  deallocate(flamewidth)
#endif
  
  return
  
end subroutine fl_flameSpeed
