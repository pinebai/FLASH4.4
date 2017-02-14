!!****if* source/physics/Gravity/GravityMain/Poisson/BHTree/Gravity_accelOneRow
!!
!! NAME
!!
!!  Gravity_accelOneRow
!!
!!
!! SYNOPSIS
!!
!!  call Gravity_accelOneRow(integer(IN)  :: pos(2),
!!                           integer(IN)  :: sweepDir,
!!                           integer(IN)  :: blockID,
!!                           integer(IN)  :: numCells,
!!                           real(INOUT)  :: grav(numCells),
!!                           integer(IN),optional :: potentialIndex,
!!                           integer(IN),optional :: extraAccelVars(MDIM))
!!
!! DESCRIPTION
!!
!!  Compute components of the zone-averaged gravitational
!!  acceleration on a mesh block.  A single vector component
!!  of the acceleration is computed.
!!
!!  This routine computes the gravitational acceleration for a row
!!  of zones in a specified direction in a given block. First-order
!!  finite-volume differencing is used everywhere.  Formulae based
!!  on long stencils (usually high-order) may produce differences
!!  at the block boundaries for siblings as hydro solver may require
!!  several valid guard cells (e.g., PPM with parabolic
!!  interpolation for force terms needs 3 valid guard cells). Not
!!  providing such valid data may result in violation of conservation. 
!!
!! ARGUMENTS
!!
!!  pos     -       Row indices transverse to the sweep direction
!!  sweepDir   -       The sweep direction:  allowed values are 
!!              SWEEP_X, SWEEP_Y, and SWEEP_Z. These values are defined
!!              in constants.h.
!!  blockID   -     The local identifier of the block to work on
!!  grav()   -       Array to receive result
!!  numCells -       Number of cells to update in grav array
!!  potentialIndex      -  if specified,  Variable # to take as potential.
!!                         Default is GPOT_VAR for the potential stored in the
!!                         gpot slot of unk, which should correspond to the
!!                         potential at the current timestep.
!!  extraAccelVars      -  if specified,  Variables from which extra accelerations
!!                         are taken. Used to identify the UNK variables
!!                         that contain sink-on-gas accelerations when
!!                         sink particles are used.
!!
!! NOTES
!!
!!  If certain variables declared by the sink particles inplementation are
!!  declared, it is assumed that sink particles are in use.
!!  The sets of variables to make this determination are
!!    o  those given by extraAccelVars   extraAccelVars if present;
!!    o  {SGXO_VAR, SGYO_VAR, SGZO_VAR}  if potentialIndex is GPOL_VAR;
!!    o  {SGAX_VAR, SGAY_VAR, SGAZ_VAR}  otherwise.
!!  If it is assumed that sink particles are in use, then the acceleration
!!  returned in the grav array will have the appropriate sink particle
!!  acceleration component added to the acceleration computed by differencing
!!  the potential variable given by potentialIndex.
!!***

!!REORDER(4): solnVec

subroutine Gravity_accelOneRow (pos, sweepDir, blockID, numCells, grav, &
                                potentialIndex, extraAccelVars)


  use Grid_interface, ONLY : Grid_getBlkPhysicalSize, Grid_getBlkPtr, &
    Grid_releaseBlkPtr, Grid_getBlkIndexLimits
  !use Driver_interface, ONLY : Driver_abortFlash
  use grv_bhInterface, ONLY : grv_accExtrnPotRow
  use Gravity_data, ONLY: grv_defaultGpotVar, grv_useExternalPotential, &
    grv_usePoissonPotential, grv_poisson_max, grv_sink_max

  implicit none

#include "Flash.h"
#ifdef Grid_releaseBlkPtr
! disabling per-block drift logging for this routine because it is called too much
#undef Grid_releaseBlkPtr
#endif

#include "constants.h"

  integer, dimension(2), intent(in) :: pos
  integer, intent(in)               :: sweepDir, blockID,  numCells
  real, intent(inout)               :: grav(numCells)
  integer, intent(IN),optional      :: potentialIndex
  integer, intent(IN),OPTIONAL      :: extraAccelVars(MDIM)

  real            :: blockSize(MDIM)
  real, POINTER, DIMENSION(:,:,:,:) :: solnVec
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer         :: ii, iimin, iimax, j, k
  real            :: gpot(numCells), delxinv, tmpdr
  real :: grav_poisson(numCells), grav_ext(numCells), grav_sink(numCells)
  integer         :: potVar, nxbBlock, nybBlock, nzbBlock
  integer         :: sink_ax_index, sink_ay_index, sink_az_index
  integer         :: pois_ax_index, pois_ay_index, pois_az_index
  !real, dimension(MDIM) :: xxx

  !==================================================
  
  
  call Grid_getBlkPhysicalSize(blockID, blockSize)

  
  call Grid_getBlkPtr(blockID, solnVec)


  call Grid_getBlkIndexLimits(blockId, blkLimits, blkLimitsGC)
  nxbBlock = blkLimits(HIGH,IAXIS) - blkLimits(LOW,IAXIS) + 1
  nybBlock = blkLimits(HIGH,JAXIS) - blkLimits(LOW,JAXIS) + 1
  nzbBlock = blkLimits(HIGH,KAXIS) - blkLimits(LOW,KAXIS) + 1

  iimin   = 1
  iimax   = numCells

  ! IF a variable index is explicitly specified, assume that as the potential
  ! otherwise use the default current potential GPOT_VAR  
  if(present(potentialIndex)) then
    potVar=potentialIndex
  else if (grv_defaultGpotVar > 0) then
    potVar=grv_defaultGpotVar
  else
    potVar=GPOT_VAR
  end if

  sink_ax_index = 0
  sink_ay_index = 0
  sink_az_index = 0

  if (present(extraAccelVars)) then
  ! select the current or the old indeces for poisson and sink accelerations
     sink_ax_index = extraAccelVars(1)
     sink_ay_index = extraAccelVars(2)
     sink_az_index = extraAccelVars(3)
  else if (potVar .eq. GPOT_VAR) then
#if defined(GACX_VAR) && defined(GACY_VAR) && defined(GACZ_VAR)
    pois_ax_index = GACX_VAR
    pois_ay_index = GACY_VAR
    pois_az_index = GACZ_VAR
#endif
#if defined(SGAX_VAR) && defined(SGAY_VAR) && defined(SGAZ_VAR)
    sink_ax_index = SGAX_VAR
    sink_ay_index = SGAY_VAR
    sink_az_index = SGAZ_VAR
#endif
  else if (potVar .eq. GPOL_VAR) then
#if defined(GAOX_VAR) && defined(GAOY_VAR) && defined(GAOZ_VAR)
    pois_ax_index = GAOX_VAR
    pois_ay_index = GAOY_VAR
    pois_az_index = GAOZ_VAR
#endif
! in some cases SG[XYZ]O are not defined; if so, let us use SGA[XYZ] for the old
! sink accelerations    
#if defined(SGAX_VAR) && defined(SGAY_VAR) && defined(SGAZ_VAR)
    sink_ax_index = SGAX_VAR
    sink_ay_index = SGAY_VAR
    sink_az_index = SGAZ_VAR
#endif
#if defined(SGXO_VAR) && defined(SGYO_VAR) && defined(SGZO_VAR)
    sink_ax_index = SGXO_VAR
    sink_ay_index = SGYO_VAR
    sink_az_index = SGZO_VAR
#endif
  endif
!  if ((potVar .ne. GPOT_VAR) .and. (potVar .ne. GPOL_VAR)) then
!     print *, "Gravity_accelOneRow called with neither GPOT_VAR nor GPOL_VAR."
!     call Driver_abortFlash("Gravity_accelOneRow called with neither GPOT_VAR nor GPOL_VAR.")
!  endif

  ! 1. Contribution to acceleration from the Poisson solver (gas)
  if (grv_usePoissonPotential) then
#ifdef GRAV_TREE_ACC
    if (sweepDir .eq. SWEEP_X) then
      grav_poisson = solnVec(pois_ax_index, :, pos(1), pos(2))
    else if (sweepDir .eq. SWEEP_Y) then
      grav_poisson = solnVec(pois_ay_index, pos(1), :, pos(2))
    else
      grav_poisson = solnVec(pois_az_index, pos(1), pos(2), :)
    endif
#else
    ! GRAV_TREE_ACC = 0; calculating potential
    ! Get row of potential values and compute inverse of zone spacing  
    if (sweepDir == SWEEP_X) then                     ! x-direction
      delxinv = real(nxbBlock) / blockSize(IAXIS)
      gpot(:) = solnVec(potVar,:,pos(1),pos(2))
    elseif (sweepDir == SWEEP_Y) then                 ! y-direction
      delxinv = real(nybBlock) / blockSize(JAXIS)
      gpot(:) = solnVec(potVar,pos(1),:,pos(2))
    else                                              ! z-direction
      delxinv = real(nzbBlock) / blockSize(KAXIS)
      gpot(:) = solnVec(potVar,pos(1),pos(2),:)
    endif
    ! first-order differences preserves conservation
    delxinv = 0.5e0 * delxinv
    do ii = iimin+1, iimax-1
      grav_poisson(ii) = delxinv * (gpot(ii-1) - gpot(ii+1))
    enddo
    grav_poisson(iimin) = grav_poisson(iimin+1)     ! this is invalid data - must not be used
    grav_poisson(iimax) = grav_poisson(iimax-1)

    ! if acceleration variables are present, copy derivative of the potential into them 
#if defined(GACX_VAR) && defined(GACY_VAR) && defined(GACZ_VAR)
    if (sweepDir == SWEEP_X) then
      solnVec(GACX_VAR, :, pos(1), pos(2)) = grav_poisson
    else if (sweepDir == SWEEP_Y) then
      solnVec(GACY_VAR, pos(1), :, pos(2)) = grav_poisson
    else if (sweepDir == SWEEP_Z) then
      solnVec(GACZ_VAR, pos(1), pos(2), :) = grav_poisson
    endif
#endif

#endif
  else
    grav_poisson = 0.0
  endif

  ! 2. Contribution to acceleration from sink particles
#if defined(SGAX_VAR) && defined(SGAY_VAR) && defined(SGAZ_VAR)
  if (sweepDir .eq. SWEEP_X) then
    grav_sink = solnVec(sink_ax_index, :, pos(1), pos(2))
  else if (sweepDir .eq. SWEEP_Y) then
    grav_sink = solnVec(sink_ay_index, pos(1), :, pos(2))
  else
    grav_sink = solnVec(sink_az_index, pos(1), pos(2), :)
  endif
#else
  grav_sink(:) = 0.0
#endif

  ! 3. Contribution to the acceleration from the external field
  call grv_accExtrnPotRow(pos, sweepDir, blockID, numCells, grav_ext)

  !! sum all contributions
  do ii = iimin, iimax
     grav(ii) = grav_poisson(ii) + grav_ext(ii) + grav_sink(ii)
     if (abs(grav_poisson(ii)) > grv_poisson_max) &
       grv_poisson_max = abs(grav_poisson(ii))
     if (abs(grav_sink(ii)) > grv_sink_max) &
       grv_sink_max = abs(grav_sink(ii))
  enddo
  !print *, "Poisson: ", grav_poisson
  !print *, "Sink: ", grav_sink
  !print *, "Extern: ", grav_ext

  call Grid_releaseBlkPtr(blockID, solnVec)
  
  return
   
end subroutine Gravity_accelOneRow


