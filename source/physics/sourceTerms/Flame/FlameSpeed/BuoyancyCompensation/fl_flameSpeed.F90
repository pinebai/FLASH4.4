!!****if* source/physics/sourceTerms/Flame/FlameSpeed/BuoyancyCompensation/fl_flameSpeed
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
!!
! Dean Townsley 2008
!

!! Several suppression options are implemented which "turn off" the flame
!! by setting the flame speed to zero, generally based on time and axial angle
!! (theta).  One acts only on the buoy comp component and one on the whole flame
!! speed.
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


subroutine fl_flameSpeed(solnData, flamespeed, blockID, nlayers)

#include "constants.h"
#include "Flash.h"
  
  use fl_fsData, ONLY :fl_fsUseConstFlameSpeed, fl_fsQuench,&
       fl_fsConstFlameSpeed, fl_fsConstFlameWidth, fl_fsUseTFI, &
       fl_fsQuenchDens0, fl_fsQuenchDens1,fl_fsQuenchInvDDens,&
       fl_fsGcdFlameSuppress, fl_fsGcdFlameSuppressTime, &
       fl_fsGcdFlameSuppressCosTheta, fl_fsBuoyCompSuppress, &
       fl_fsBuoyCompSuppressTime, fl_fsBuoyCompSuppressCosTheta, &
       fl_fsGeom, fl_fsMDelta
  use Driver_interface, ONLY : Driver_abortFlash, Driver_getSimTime
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getCellCoords, &
       Grid_getDeltas
  use fl_fsAtwoodInterface, ONLY : fl_fsAtwood
  use fl_fsLaminarInterface, ONLY : fl_fsLaminarFlameSpeedBlock
  use fl_fsTFIInterface, ONLY : fl_fsTFIFlameSpeedBlock
  use Gravity_interface, only : Gravity_accelOneBlock

  implicit none
  
  integer, intent(in) :: blockID
  real,dimension(:,:,:,:),pointer :: solnData
  real,dimension(:,:,:),intent(out) :: flamespeed
  integer, intent(in) :: nlayers

  integer :: sizeI, sizeJ, sizeK
  integer :: istat
  
  integer, dimension(LOW:HIGH,MDIM)     :: blkLimits, blkLimitsGC, compLimits
  real, dimension(MDIM)                 :: deltas

  real, allocatable, dimension(:,:,:)   :: atwood,dens,flamewidth,quench_limit
  real, dimension(:,:,:,:), allocatable :: gvec
  real,allocatable,dimension(:) :: iCoord,jCoord,kCoord

  integer :: i,j,k
  real :: time, costheta, bcspeed, magg, x, dx

  ! allow jumpout with constant flame speed
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

  sizeI=blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  sizeJ=blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
  sizeK=blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

  ! make an array of limits over which computation is being performed
  compLimits(LOW,IAXIS) = blkLimits(LOW,IAXIS)-nlayers
  compLimits(HIGH,IAXIS) = blkLimits(HIGH,IAXIS)+nlayers
  compLimits(LOW,JAXIS) = blkLimits(LOW,JAXIS)-K2D*nlayers
  compLimits(HIGH,JAXIS) = blkLimits(HIGH,JAXIS)+K2D*nlayers
  compLimits(LOW,KAXIS) = blkLimits(LOW,KAXIS)-K3D*nlayers
  compLimits(HIGH,KAXIS) = blkLimits(HIGH,KAXIS)+K3D*nlayers

  allocate(flamewidth(sizeI,sizeJ,sizeK),STAT=istat)
  if (istat /= 0) call Driver_abortFlash("Cannot allocate flamewidth in fl_flameSpeed")

  ! allow jumpout with constant flame speed with TFI
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
     deallocate(flamewidth)
     return
  endif
 
  call Driver_getSimTime(time)

  allocate(atwood(sizeI,sizeJ,sizeK),STAT=istat)
  if (istat /= 0) call Driver_abortFlash("Cannot allocate atwood in fl_flameSpeed")
  allocate(dens(sizeI,sizeJ,sizeK),STAT=istat)
  if (istat /= 0) call Driver_abortFlash("Cannot allocate dens in fl_flameSpeed")
  allocate(gvec(NDIM,sizeI,sizeJ,sizeK), stat=istat)
  if (istat/=0) call Driver_abortFlash("Failed to allocate gvec in fl_flameSpeed")
  if (fl_fsUseTFI) then
     allocate(quench_limit(sizeI,sizeJ,sizeK),STAT=istat)
     if (istat /= 0) call Driver_abortFlash("Cannot allocate quench_limit in fl_flameSpeed")
  endif
  
  ! get necessary coordinates depending on options and geometry
  if (fl_fsBuoyCompSuppress .or. fl_fsGcdFlameSuppress) then
     select case (fl_fsGeom)
     case (CARTESIAN)
        allocate(iCoord(sizeI),STAT=istat)
        if (istat /= 0) call Driver_abortFlash("Cannot allocate iCoord in fl_flameSpeed")
        allocate(jCoord(sizeJ),STAT=istat)
        if (istat /= 0) call Driver_abortFlash("Cannot allocate jCoord in fl_flameSpeed")
        allocate(kCoord(sizeK),STAT=istat)
        if (istat /= 0) call Driver_abortFlash("Cannot allocate kCoord in fl_flameSpeed")
        call Grid_getCellCoords(IAXIS,blockID,CENTER,.true.,iCoord,sizeI)
        call Grid_getCellCoords(JAXIS,blockID,CENTER,.true.,jCoord,sizeJ)
        call Grid_getCellCoords(KAXIS,blockID,CENTER,.true.,kCoord,sizeK)
     case (CYLINDRICAL)
        allocate(iCoord(sizeI),STAT=istat)
        if (istat /= 0) call Driver_abortFlash("Cannot allocate iCoord in fl_flameSpeed")
        allocate(jCoord(sizeJ),STAT=istat)
        if (istat /= 0) call Driver_abortFlash("Cannot allocate jCoord in fl_flameSpeed")
        call Grid_getCellCoords(IAXIS,blockID,CENTER,.true.,iCoord,sizeI)
        call Grid_getCellCoords(JAXIS,blockID,CENTER,.true.,jCoord,sizeJ)
     case default
        call Driver_abortFlash("bad geometry for angle-based suppression")
     end select
  endif

  ! this returns laminar flame speed and estimate of unburned density at this pressure
  call fl_fsLaminarFlameSpeedBlock(solnData, flamespeed, dens,  &
                                   compLimits, flamewidth)
  ! this accounts for Turbulence-Flame Interaction
  if (fl_fsUseTFI) &
    call fl_fsTFIFlameSpeedBlock(solnData, flamespeed, flamewidth, dx, &
      compLimits, quench_limit)
  ! get atwood number estimate
  call fl_fsAtwood(solnData, atwood, dens, compLimits)
  ! calculate gravitational acceleration
  call Gravity_accelOneBlock(blockID, nlayers, gvec)

  ! now calculate buoyancy-compensated flame speed
  do k = compLimits(LOW,KAXIS), compLimits(HIGH,KAXIS)
     do j = compLimits(LOW,JAXIS), compLimits(HIGH,JAXIS)
        do i = compLimits(LOW,IAXIS), compLimits(HIGH,IAXIS)

           ! calculate magnitude of g
           magg = gvec(IAXIS,i,j,k)**2
           if (NDIM>=2) magg = magg+gvec(JAXIS,i,j,k)**2
           if (NDIM==3) magg = magg+gvec(KAXIS,i,j,k)**2
           magg = sqrt(magg)

           ! for these options we need costheta
           if (fl_fsBuoyCompSuppress .or. fl_fsGcdFlameSuppress) then
              select case (fl_fsGeom)
              case (CARTESIAN)
                costheta = kCoord(k)/sqrt(iCoord(i)**2+jCoord(j)**2+kCoord(k)**2)
              case (CYLINDRICAL)
                 costheta = jCoord(j)/sqrt(iCoord(i)**2+jCoord(k)**2)
              case default
                 call Driver_abortFlash("bad geometry for angle-based suppression")
              end select
           endif

           ! this is the 0.5*sqrt(A g m Delta) flame speed
           if (fl_fsBuoyCompSuppress .and.  &
               time > fl_fsBuoyCompSuppressTime .and. costheta <= fl_fsBuoyCompSuppressCosTheta ) then
              bcspeed = 0.0
           else
              bcspeed = 0.5*sqrt( atwood(i,j,k)*magg*fl_fsMDelta)
           endif

           if (fl_fsUseTFI) then
              if (quench_limit(i,j,k) > 0.0) then
                 bcspeed = min( quench_limit(i,j,k), bcspeed)
              endif
           endif

           ! recall flame speed is currently set to laminar (or turbluent) speed
           flamespeed(i,j,k) = max( flamespeed(i,j,k), bcspeed)

           ! quench flame at low densities if requested
           if (fl_fsQuench) then
              if (dens(i,j,k) < fl_fsQuenchDens0) then
                 flamespeed(i,j,k) = 0.0
              else if (dens(i,j,k) < fl_fsQuenchDens1) then
                 x = ( dens(i,j,k) -fl_fsQuenchDens0)*fl_fsQuenchInvDDens
                 flamespeed(i,j,k) = (3.0*x-2.0*x**2)*x*flamespeed(i,j,k)
              endif
           endif

           !  Suppress flame in selected region
           if ( fl_fsGcdFlameSuppress .and. &
                time > fl_fsGcdFlameSuppressTime .and. costheta <= fl_fsGcdFlameSuppressCosTheta ) then
                flamespeed(i,j,k) = 0.0
           endif
          
        enddo
     enddo
  enddo

  ! save flame speed if variable is provided
#ifdef FSPD_VAR
  solndata(FSPD_VAR, compLimits(LOW,IAXIS):compLimits(HIGH,IAXIS), &
                     compLimits(LOW,JAXIS):compLimits(HIGH,JAXIS), &
                     compLimits(LOW,KAXIS):compLimits(HIGH,KAXIS)) &
                      = flamespeed(compLimits(LOW,IAXIS):compLimits(HIGH,IAXIS), &
                                   compLimits(LOW,JAXIS):compLimits(HIGH,JAXIS), &
                                   compLimits(LOW,KAXIS):compLimits(HIGH,KAXIS))
#endif
  
  deallocate(atwood)
  deallocate(dens)
  deallocate(flamewidth)
  deallocate(gvec)
  if (fl_fsUseTFI) then
     deallocate(quench_limit)
  endif

  if (fl_fsBuoyCompSuppress .or. fl_fsGcdFlameSuppress) then
     select case (fl_fsGeom)
     case (CARTESIAN)
        deallocate(iCoord)
        deallocate(jCoord)
        deallocate(kCoord)
     case (CYLINDRICAL)
        deallocate(iCoord)
        deallocate(jCoord)
     end select
  endif

  return
  
end subroutine fl_flameSpeed
