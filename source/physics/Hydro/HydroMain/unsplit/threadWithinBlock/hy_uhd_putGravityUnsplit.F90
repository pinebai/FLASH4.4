!!****if* source/physics/Hydro/HydroMain/unsplit/threadWithinBlock/hy_uhd_putGravityUnsplit
!!
!! NAME
!!
!!  hy_uhd_PutGravityUnsplit
!!
!! SYNOPSIS
!!
!!  hy_uhd_putGravityUnsplit( integer(IN) :: blockID,
!!                            integer(IN) :: blkLimitsGC(2,MDIM),
!!                            integer(IN) :: dataSize(MDIM),
!!                            real   (IN) :: dt,
!!                            real   (IN) :: dtOld,
!!                            real(OUT)   :: gravX(:,:,:),
!!                            real(OUT)   :: gravY(:,:,:),
!!                            real(OUT)   :: gravZ(:,:,:))
!!
!! ARGUMENTS
!!
!!  blockID     - a current block ID
!!  blkLimitsGC - an array that holds the lower and upper indices of the section
!!                of block with the guard cells
!!  dataSize    - dimensions for gravX, gravY and gravZ arrays
!!  dt          - timestep
!!  dtOld       - old timestep (needed for temporal extrapolations of gravity)
!!  gravX       - gravity components in x-direcition at time steps n
!!  gravY       - gravity components in y-direcition at time steps n
!!  gravZ       - gravity components in z-direcition at time steps n
!!
!! DESCRIPTION
!!
!!  This routine puts gravity components to arrays gravX, gravY and gravZ:
!!  gravX(:,:,:) includes gravity components at time step n
!!
!!*** 

Subroutine hy_uhd_putGravityUnsplit&
     (blockID,blkLimitsGC,dataSize,dt,dtOld,gravX,gravY,gravZ, potentialIndex, lastCall)

  use Gravity_interface, ONLY : Gravity_accelOneRow
  use Hydro_data, ONLY : hy_threadWithinBlock

  use Hydro_data, ONLY: hy_gpotVar, hy_extraAccelVars

  implicit none

#include "Flash.h"
#include "constants.h"
#include "UHD.h"
#include "Flash_omp.h"

  !! ---- Argument List ----------------------------------
  integer, intent(IN) :: blockID
  integer, dimension(LOW:HIGH,MDIM), intent(IN) :: blkLimitsGC
  integer, dimension(MDIM), intent(IN) :: dataSize
  real,    intent(IN) :: dt, dtOld

#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC), intent(OUT) :: &
       gravX,gravY,gravZ
#else
  real, dimension(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)), intent(OUT) :: &
       gravX,gravY,gravZ
#endif
  integer, intent(IN), OPTIONAL :: potentialIndex
  logical, intent(IN), OPTIONAL :: lastCall
  !! -----------------------------------------------------

  integer, dimension(2) :: gravPos
  integer :: ix, iy, iz
  integer :: potVar
  logical :: tailCall

  if (present(lastCall)) then
     tailCall = lastCall
  else
     tailCall = .FALSE.
  end if

  !$omp parallel if (hy_threadWithinBlock) &
  !$omp default(none) &
  !$omp shared(potentialIndex,tailCall,hy_gpotVar,hy_extraAccelVars,&
  !$omp gravX,gravY,gravZ,blockID,dataSize,blkLimitsGC) &
  !$omp private(ix,iy,iz,gravPos,potVar)


#if NDIM == 3
  !$omp do schedule(static)
#endif
  do iz=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
#if NDIM == 2
     !$omp do schedule(static)
#endif
     do iy=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        gravPos(1)=iy
        gravPos(2)=iz
        if (tailCall .AND. hy_gpotVar > 0) then
           if (present(potentialIndex)) then
              potVar = potentialIndex
           else
              potVar = hy_gpotVar
           end if
           if (hy_extraAccelVars(1)>0) then
              call Gravity_accelOneRow(gravPos,DIR_X,blockID,dataSize(IAXIS),gravX(:,iy,iz),potVar,&
                                       extraAccelVars=hy_extraAccelVars)
           else
              call Gravity_accelOneRow(gravPos,DIR_X,blockID,dataSize(IAXIS),gravX(:,iy,iz),potVar)
           end if
        else if (present(potentialIndex)) then
           call Gravity_accelOneRow(gravPos,DIR_X,blockID,dataSize(IAXIS),gravX(:,iy,iz),potentialIndex)
        else
#if defined(GPOT_VAR) && defined(FLASH_GRAVITY_TIMEDEP)
        ! Gravity implementation defines FLASH_GRAVITY_TIMEDEP -> time-dependent gravity field
        ! gravity at time step n
           call Gravity_accelOneRow(gravPos,DIR_X,blockID,dataSize(IAXIS),gravX(:,iy,iz),GPOT_VAR)
#else
        ! FLASH_GRAVITY_TIMEDEP not defined -> assume time-independent gravity field.
        ! Also if GPOT_VAR is defined -> use current accel without time
        ! interpolation, i.e., handle like time-independent gravity field - KW
           call Gravity_accelOneRow(gravPos,DIR_X,blockID,dataSize(IAXIS),gravX(:,iy,iz))
#endif
        endif
     enddo
#if NDIM == 2
     !$omp end do
#endif
  enddo
#if NDIM == 3
  !$omp end do
#endif


  if (NDIM >= 2) then
#if NDIM == 3
     !$omp do schedule(static)
#endif
     do iz=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
#if NDIM == 2
        !$omp do schedule(static)
#endif
        do ix=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
           gravPos(1)=ix
           gravPos(2)=iz
           if (tailCall .AND. hy_gpotVar > 0) then
              if (present(potentialIndex)) then
                 potVar = potentialIndex
              else
                 potVar = hy_gpotVar
              end if
              if (hy_extraAccelVars(2)>0) then
                 call Gravity_accelOneRow(gravPos,DIR_Y,blockID,dataSize(JAXIS),gravY(ix,:,iz),potVar,&
                                          extraAccelVars=hy_extraAccelVars)
              else
                 call Gravity_accelOneRow(gravPos,DIR_Y,blockID,dataSize(JAXIS),gravY(ix,:,iz),potVar)
              end if
           else if (present(potentialIndex)) then
              call Gravity_accelOneRow(gravPos,DIR_Y,blockID,dataSize(JAXIS),gravY(ix,:,iz),potentialIndex)
           else
#if defined(GPOT_VAR) && defined(FLASH_GRAVITY_TIMEDEP)
           ! gravity at time step n
              call Gravity_accelOneRow(gravPos,DIR_Y,blockID,dataSize(JAXIS),gravY(ix,:,iz),GPOT_VAR)
#else
              call Gravity_accelOneRow(gravPos,DIR_Y,blockID,dataSize(JAXIS),gravY(ix,:,iz))
#endif
           end if
        enddo
#if NDIM == 2
     !$omp end do
#endif
     enddo
#if NDIM == 3
  !$omp end do
#endif


     if (NDIM == 3) then
        !$omp do schedule(static)
        do iy=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
           do ix=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
              gravPos(1)=ix
              gravPos(2)=iy
              if (tailCall .AND. hy_gpotVar > 0) then
                 if (present(potentialIndex)) then
                    potVar = potentialIndex
                 else
                    potVar = hy_gpotVar
                 end if
                 if (hy_extraAccelVars(3)>0) then
                    call Gravity_accelOneRow(gravPos,DIR_Z,blockID,dataSize(KAXIS),gravZ(ix,iy,:),potVar,&
                                             extraAccelVars=hy_extraAccelVars)
                 else
                    call Gravity_accelOneRow(gravPos,DIR_Z,blockID,dataSize(KAXIS),gravZ(ix,iy,:),potVar)
                 end if
              else if (present(potentialIndex)) then
                 call Gravity_accelOneRow(gravPos,DIR_Z,blockID,dataSize(KAXIS),gravZ(ix,iy,:),potentialIndex)
              else
#if defined(GPOT_VAR) && defined(FLASH_GRAVITY_TIMEDEP)
              ! gravity at time step n
                 call Gravity_accelOneRow(gravPos,DIR_Z,blockID,dataSize(KAXIS),gravZ(ix,iy,:),GPOT_VAR)
#else
                 call Gravity_accelOneRow(gravPos,DIR_Z,blockID,dataSize(KAXIS),gravZ(ix,iy,:))
#endif
              end if
           enddo
        enddo
        !$omp end do
     endif
  endif
  !$omp end parallel

End Subroutine hy_uhd_putGravityUnsplit
