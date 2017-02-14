!!****if* source/physics/Hydro/HydroMain/unsplit_rad/threadWithinBlock/hy_uhd_getFaceFlux
!!
!! NAME
!!
!!  hy_uhd_getFaceFlux
!!
!! SYNOPSIS
!!
!!  hy_uhd_getFaceFlux( integer(IN) :: blockID,
!!                      integer(IN) :: blkLimits(2,MDIM),
!!                      integer(IN) :: blkLimitsGC(2,MDIM), 
!!                      integer(IN) :: datasize(MDIM),
!!                      real(IN)    :: del(MDIM),
!!                      real(OUT)   :: xflux(:,:,:,:),
!!                      real(OUT)   :: yflux(:,:,:,:),
!!                      real(OUT)   :: zflux(:,:,:,:),
!!                      real, pointer, dimension(:,:,:,:) :: scrchFaceXPtr,
!!                      real, pointer, dimension(:,:,:,:) :: scrchFaceYPtr,
!!                      real, pointer, dimension(:,:,:,:) :: scrchFaceZPtr,
!!                      real, pointer, dimension(:,:,:,:) :: scrch_Ptr,
!!                      real, pointer, optional, dimension(:,:,:,:,:) :: hy_SpcR,hy_SpcL,hy_SpcSig,
!!                      logical,optional(IN) :: lastCall )
!!
!! ARGUMENTS
!!
!!  blockID           - a current block ID
!!  blkLimits         - an array that holds the lower and upper indices of the section
!!                      of block without the guard cells
!!  blkLimitsGC       - an array that holds the lower and upper indices of the section
!!                      of block with the guard cells
!!  datasize          - data size for boundary extrapolated data, boundary data
!!  del               - grid deltas
!!  xflux,yflux,zflux - face fluxes at each {x,y,z} direction
!!  scrchFaceXPtr,scrchFaceYPtr,scrchFaceZPtr - Pointers to the scrch array (for left/right states and faces)
!!  scrch_Ptr         - Pointer to the scrch array (SCRATCH_CTR)
!!  hy_SpcR,hy_SpcL,hy_SpcSig - Pointers for Species and mass scalar recon.
!!  lastCall          - if true then store flux data in scratch array
!!
!! DESCRIPTION
!!
!!  This routine computes high-order Godunov fluxes at cell interface centers 
!!  for each spatial direction using a choice of Riemann solvers.
!!  Choices of Riemann solvers are Roe-type, HLL(E), HLLC, HLLD, 
!!  Hybrid (HLLC+HLL for hydro; HLLD+HLL for MHD), local Lax-Friedrichs, and Marquina solvers.
!!
!!*** 

!!REORDER(4):U, scrchFaceXPtr,scrchFaceYPtr,scrchFaceZPtr, scrch_Ptr, [xyz]flux

!#define COMPUTE_DT_FLUX

#include "Flash.h"
subroutine hy_uhd_getFaceFlux ( blockID,blkLimits,blkLimitsGC,datasize,del,&
                                xflux,yflux,zflux,&
                                scrchFaceXPtr,scrchFaceYPtr,scrchFaceZPtr,scrch_Ptr,&
                                hy_SpcR,hy_SpcL,hy_SpcSig,lastCall)

  use Hydro_data,                    ONLY : hy_nref,hy_kref,          &
                                            hy_RiemannSolver,         &
                                            hy_useDiffuse,            &
                                            hy_useViscosity,          &
                                            hy_useConductivity,       &
                                            hy_use_avisc,             &
                                            hy_cvisc,                 &
                                            hy_updateHydroFluxes,     &
                                            hy_addThermalFlux,        &
                                            hy_geometry,              &
                                            hy_numXN,                 &
                                            hy_fullRiemannStateArrays,         &
                                            hy_fullSpecMsFluxHandling,&
                                            hy_useAuxEintEqn,         &
                                            hy_hydroComputeDtOption,  &
                                            hy_threadWithinBlock,     &
                                            hy_EOSforRiemann

  use hy_uhd_interface,              ONLY : hy_uhd_addViscousFluxes,  &
                                            hy_uhd_addThermalFluxes,  &
                                            hy_uhd_Roe, &
                                            hy_uhd_LLF, &
                                            hy_uhd_HLL, &
                                            hy_uhd_HLLC,&
                                            hy_uhd_Marquina,&
                                            hy_uhd_MarquinaModified,&
                                            hy_uhd_setMinTimeStep

  use Grid_interface,                ONLY : Grid_getBlkPtr, &
                                            Grid_releaseBlkPtr, &
                                            Grid_getCellCoords
  use Conductivity_interface,        ONLY : Conductivity
  use Viscosity_interface,           ONLY : Viscosity
  use Timers_interface,              ONLY : Timers_start, Timers_stop

#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
  use Hydro_data,                    ONLY : hy_mref, &
                                            hy_forceHydroLimit,&
                                            hy_useResistivity, &
                                            hy_useMagneticResistivity,&
                                            hy_useBiermann, &
                                            hy_biermannSource, &
                                            hy_useBiermann1T
  use hy_uhd_interface,              ONLY : hy_uhd_addResistiveFluxes, &
                                            hy_uhd_HLLD,&
                                            hy_uhd_addBiermannBatteryTerms
  use MagneticResistivity_interface, ONLY : MagneticResistivity
#endif


  implicit none

#include "constants.h"
#include "Eos.h"
#include "UHD.h"
#include "Flash_omp.h"

  !! Arguments type declaration ------------------------------
  integer, intent(IN)  :: blockID
  integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimits, blkLimitsGC
  integer, dimension(MDIM), intent(IN)         :: datasize
  real,    dimension(MDIM), intent(IN)         :: del

#ifdef FIXEDBLOCKSIZE
  real, intent(OUT) :: xflux (NFLUXES,&
       GRID_ILO_GC:GRID_IHI_GC, &
       GRID_JLO_GC:GRID_JHI_GC, &
       GRID_KLO_GC:GRID_KHI_GC)
  real, intent(OUT) :: yflux (NFLUXES,&
       GRID_ILO_GC:GRID_IHI_GC, &
       GRID_JLO_GC:GRID_JHI_GC, &
       GRID_KLO_GC:GRID_KHI_GC)
  real, intent(OUT) :: zflux (NFLUXES,&
       GRID_ILO_GC:GRID_IHI_GC, &
       GRID_JLO_GC:GRID_JHI_GC, &
       GRID_KLO_GC:GRID_KHI_GC)
#else
  real, intent(OUT) :: xflux (NFLUXES,&
      blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
      blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
      blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))
  real, intent(OUT) :: yflux (NFLUXES,&
      blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
      blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
      blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))
  real, intent(OUT) :: zflux (NFLUXES,&
      blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
      blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
      blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))
#endif
  real, pointer, dimension(:,:,:,:) :: scrchFaceXPtr,scrchFaceYPtr,scrchFaceZPtr
  real, pointer, dimension(:,:,:,:) :: scrch_Ptr
  real, pointer, optional, dimension(:,:,:,:,:) :: hy_SpcR,hy_SpcL,hy_SpcSig
  logical, optional, intent(IN) :: lastCall
  !! ---------------------------------------------------------

  integer :: i0,imax,j0,jmax,k0,kmax,jbeg,jend,kbeg,kend, i,j,k, ierr
  real, pointer, dimension(:,:,:,:) :: U
  !U contains (DENS,VELX,VELY,VELZ,(MAGX,MAGY,MAGZ),PRES + GAMC,GAME,EINT,TEMP)
  real, dimension(HY_VARINUMMAX) :: VL,VR 
  real, dimension(NSPECIES)    :: speciesArr

#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_ILO_GC:GRID_IHI_GC, &
                  GRID_JLO_GC:GRID_JHI_GC, &
                  GRID_KLO_GC:GRID_KHI_GC) &
                  :: viscDynamic,cond

#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
  real, dimension(GRID_ILO_GC:GRID_IHI_GC, &
                  GRID_JLO_GC:GRID_JHI_GC, &
                  GRID_KLO_GC:GRID_KHI_GC) &
                  :: magResist
#endif

#else
  real, dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
                  blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
                  blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) &
                  :: viscDynamic,cond

#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
  real, dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
                  blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
                  blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) &
                  :: magResist
#endif

#endif
  real    :: viscKinematicUnused,dcffUnused
  real    :: cvisc
  integer :: k2,k3,kGrav,kUSM
  real    :: presPlus, presMinus
  real    :: weightPlus, weightMinus, weightSum

#ifdef FLASH_UHD_HYDRO
  !CD: Add these variables so we can maintain a single omp parallel
  !directive for unsplit hydro and MHD simulations.  They are not used.
  real, dimension(1,1,1) :: magResist
  real, save :: hy_mref = 0
  logical, save :: hy_useResistivity = .false.
  logical, save :: hy_useMagneticResistivity = .false.
  logical, save :: hy_useBiermann = .false.
  logical, save :: hy_useBiermann1T = .false.
  logical, save :: hy_biermannSource = .false.
#endif

#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_IHI_GC) :: xCenter, xLeft, xRight
  real, dimension(GRID_JHI_GC) :: yCenter
#else
  real, dimension(datasize(IAXIS)) :: xCenter, xLeft, xRight
  real, dimension(datasize(JAXIS)) :: yCenter
#endif
  real :: dPdr, rvol, alpha
  real, allocatable :: xcent(:), ycent(:), zcent(:)
  real :: speed, dy, dz
  integer :: ispu,isph
  integer, dimension(MDIM) :: normDir,tranDir

  real :: t_start

  xflux = 0.
  if (NDIM < 3) then
     zflux = 0.
     if (NDIM < 2) yflux = 0.
  endif

  kGrav=0
#ifdef GRAVITY
  kGrav=1
#endif

  kUSM = 0
#ifdef FLASH_USM_MHD
  if (.NOT. hy_forceHydroLimit) kUSM = 1
#endif

  !! initialize

  call Grid_getBlkPtr(blockID,U,CENTER)

  i0   = blkLimits(LOW, IAXIS)
  imax = blkLimits(HIGH,IAXIS)
  j0   = blkLimits(LOW, JAXIS)
  jmax = blkLimits(HIGH,JAXIS)
  k0   = blkLimits(LOW, KAXIS)
  kmax = blkLimits(HIGH,KAXIS)

  if (NDIM == 1) then
     jbeg= 3
     jend=-1
     kbeg= 3
     kend=-1
     k2 = 0
     k3 = 0
  elseif (NDIM == 2) then
     jbeg= j0
     jend= jmax
     kbeg= 3
     kend=-1
     k2 = 1
     k3 = 0
  else
     jbeg= j0
     jend= jmax
     kbeg= k0
     kend= kmax
     k2 = 1
     k3 = 1
  endif

  !$omp parallel if (hy_threadWithinBlock) &
  !$omp default(none) &
  !$omp shared(hy_useDiffuse,hy_mref,hy_nref,hy_useViscosity,hy_useConductivity,&
  !$omp hy_useMagneticResistivity,hy_updateHydroFluxes,&
  !$omp hy_cvisc,hy_RiemannSolver,hy_use_avisc,blkLimits,blkLimitsGC,del,&
  !$omp blockID,xflux,yflux,zflux,hy_EOSforRiemann,&
  !$omp magResist,viscDynamic,cond,lastCall,&
  !$omp i0,imax,j0,jmax,k0,kmax,jbeg,jend,kbeg,kend,U,&
  !$omp scrchFaceXPtr,scrchFaceYPtr,scrchFaceZPtr,scrch_Ptr,&
  !$omp hy_useBiermann,hy_useBiermann1T,hy_biermannSource,hy_addThermalFlux,&
  !$omp kUSM,k3,k2,hy_spcl,hy_spcr,hy_hydrocomputedtoption,hy_geometry,dataSize,kGrav,&
  !$omp hy_fullRiemannStateArrays,hy_fullSpecMsFluxHandling, &
  !$omp hy_useauxeinteqn,xCenter,yCenter,xLeft,xRight,normDir,tranDir) &
  !$omp private(i,j,k,speciesArr,VL,VR,cvisc,t_start,&
  !$omp viscKinematicUnused,dcffUnused,&
  !$omp weightPlus,weightMinus,weightSum,&
  !$omp dPdr,rvol,alpha,presPlus,presMinus,ierr,isph,speed,dy,dz)

  
  
  if (hy_useDiffuse) then
     ! call Timers_start("get diffusion")
     !! Initialize
!!$#ifdef FLASH_USM_MHD
!!$     magResist     = 0.0
!!$#endif
!!$     viscDynamic   = 0.0
!!$     cond          = 0.0

#if NDIM == 3
     !$omp do schedule(static)
#endif
     do k=kbeg-2,kend+2
#if NDIM == 2
        !$omp do schedule(static)
#endif
        do j=jbeg-2,jend+2
#if NDIM == 1
           !$omp do schedule(static)
#endif
           do i=i0-2,imax+2

              !! copy species to a temporary array
              speciesArr(:) = U(SPECIES_BEGIN:SPECIES_END,i,j,k)

              viscDynamic(i,j,k) = 0.0

              if (hy_useViscosity) then
                 !! Get viscosity
                 call Viscosity(U(TEMP_VAR,i,j,k),U(DENS_VAR,i,j,k),&
                      speciesArr,viscDynamic(i,j,k),viscKinematicUnused)
              endif

              if (hy_useConductivity .and. hy_addThermalFlux) then
                 !! Get heat conductivity
                 call Conductivity(U(TEMP_VAR,i,j,k),U(DENS_VAR,i,j,k),&
                      speciesArr,cond(i,j,k),dcffUnused,2)
              endif
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
              if (hy_useMagneticResistivity) then
                 !! Get magnetic viscosity
                 call MagneticResistivity(U(:,i,j,k),magResist(i,j,k))
              endif
#endif
           enddo
#if NDIM == 1
           !$omp end do
#endif
        enddo
#if NDIM == 2
        !$omp end do
#endif
     enddo
#if NDIM == 3
     !$omp end do
#endif

     !! For non-ideal fluxes
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
     magResist  = magResist/hy_mref
#endif
     viscDynamic= viscDynamic/hy_nref
  endif

  if (hy_geometry /= CARTESIAN) then
     !$omp single
     call Grid_getCellCoords(IAXIS,blockID, CENTER,    .true.,xCenter, dataSize(IAXIS))
     call Grid_getCellCoords(JAXIS,blockID, CENTER,    .true.,yCenter, dataSize(JAXIS))
     call Grid_getCellCoords(IAXIS,blockID, LEFT_EDGE, .true.,xLeft,   dataSize(IAXIS))
     call Grid_getCellCoords(IAXIS,blockID, RIGHT_EDGE,.true.,xRight,  dataSize(IAXIS))
     !$omp end single
  endif

 !$omp single
  normDir=0; normDir(DIR_X)=1
  tranDir=2; tranDir(DIR_X)=0
  tranDir=tranDir*kUSM
  tranDir(2)=tranDir(2)*k2
  tranDir(3)=tranDir(3)*k3
  !$omp end single
  
  !! Compute intercell fluxes using the updated left & right states
  !! Calculate x-flux first
#if NDIM == 3
  !$omp do schedule(static)
#endif  
  do k=blkLimits(LOW,KAXIS)-tranDir(DIR_Z),blkLimits(HIGH,KAXIS)+tranDir(DIR_Z)+normDir(DIR_Z)
#if NDIM == 2
  !$omp do schedule(static)
#endif  
     do j=blkLimits(LOW,JAXIS)-tranDir(DIR_Y),blkLimits(HIGH,JAXIS)+tranDir(DIR_Y)+normDir(DIR_Y)
#if NDIM == 1
  !$omp do schedule(static)
#endif  
        do i=blkLimits(LOW,IAXIS)-tranDir(DIR_X),blkLimits(HIGH,IAXIS)+tranDir(DIR_X)+normDir(DIR_X)
        
           if (hy_updateHydroFluxes) then
              VL(HY_DENS:HY_END_VARS-kGrav)=scrchFaceXPtr(HY_P01_FACEXPTR_VAR:HY_P01_FACEXPTR_VAR+HY_SCRATCH_NUM-1,i-1,j,k)
              VR(HY_DENS:HY_END_VARS-kGrav)=scrchFaceXPtr(HY_N01_FACEXPTR_VAR:HY_N01_FACEXPTR_VAR+HY_SCRATCH_NUM-1,i,  j,k)
 
#ifdef BDRY_VAR
              ! solid internal boundary
              ! Cell i and i-1:
              if (U(BDRY_VAR,i,j,k) > 0.0 .and. U(BDRY_VAR,i-1,j,k) < 0.0) then
                 VR(HY_DENS:HY_END_VARS-kGrav) = VL(HY_DENS:HY_END_VARS-kGrav)
                 VR(HY_VELX) = -VL(HY_VELX)
              end if

              if (U(BDRY_VAR,i,j,k) < 0.0 .and. U(BDRY_VAR,i-1,j,k) > 0.0) then
                 VL(HY_DENS:HY_END_VARS-kGrav) = VR(HY_DENS:HY_END_VARS-kGrav)
                 VL(HY_VELX) = -VR(HY_VELX)
              end if
#endif

              if (hy_EOSforRiemann) then
                 ! Construct initial guess for temperature
                 t_start = 0.5*(U(TEMP_VAR,i-1,j,k) + U(TEMP_VAR,i,j,k))
                 if (hy_fullSpecMsFluxHandling .AND. hy_numXN > 0) then
                    call faceStatesEos(VL,t_start,hy_SpcL(:,i  ,j,k,DIR_X))
                    call faceStatesEos(VR,t_start,hy_SpcR(:,i-1,j,k,DIR_X))
                 else
                    call faceStatesEos(VL,t_start)
                    call faceStatesEos(VR,t_start)
                 end if
              end if

              if (hy_RiemannSolver == ROE) then
                 call hy_uhd_Roe(DIR_X,VL,VR,xflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM,i,j,k),speed,ierr)
                 if(ierr /= 0) call do_error("hy_uhd_Roe", VL, VR, i,j,k, DIR_X)

              elseif (hy_RiemannSolver == HLL) then
                 call hy_uhd_HLL(DIR_X,VL,VR,xflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM,i,j,k),speed,ierr)
                 if(ierr /= 0) call do_error("hy_uhd_HLL", VL, VR, i,j,k, DIR_X)

              elseif (hy_RiemannSolver == HLLC) then
                 call hy_uhd_HLLC(DIR_X,VL,VR,xflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM,i,j,k),speed,ierr)
                 if(ierr /= 0) call do_error("hy_uhd_HLLC", VL, VR, i,j,k, DIR_X)

#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
              elseif (hy_RiemannSolver == HLLD) then
                 call hy_uhd_HLLD(DIR_X,VL,VR,xflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM,i,j,k),speed,ierr)
                 if(ierr /= 0) call do_error("hy_uhd_HLLD", VL, VR, i,j,k, DIR_X)
#endif
              elseif (hy_RiemannSolver == MARQ) then
                 call hy_uhd_Marquina(DIR_X,VL,VR,xflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM,i,j,k),speed,ierr)
                 if(ierr /= 0) call do_error("hy_uhd_Marquina", VL, VR, i,j,k, DIR_X)

              elseif (hy_RiemannSolver == MARM) then
                 call hy_uhd_MarquinaModified(DIR_X,VL,VR,xflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM,i,j,k),speed,ierr)
                 if(ierr /= 0) call do_error("hy_uhd_Marquina", VL, VR, i,j,k, DIR_X)

              elseif (hy_RiemannSolver == LLF) then
                 call hy_uhd_LLF(DIR_X,VL,VR,xflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM,i,j,k),speed,ierr)
                 if(ierr /= 0) call do_error("hy_uhd_LLF", VL, VR, i,j,k, DIR_X)

#ifdef SHOK_VAR
              elseif (hy_RiemannSolver == HYBR) then
                 if (U(SHOK_VAR,i-1,j,k) + U(SHOK_VAR,i,j,k) > 0.) then
                    !! use diffusive HLL solver for a local strong shock/rarefaction region
                    call hy_uhd_HLL(DIR_X,VL,VR,xflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM,i,j,k),speed,ierr)
                    if(ierr /= 0) call do_error("hy_uhd_HLL", VL, VR, i,j,k, DIR_X)
                 else
#ifdef FLASH_UHD_HYDRO
                    ! Hydro
                    call hy_uhd_HLLC(DIR_X,VL,VR,xflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM,i,j,k),speed,ierr)
                    if(ierr /= 0) call do_error("hy_uhd_HLLC", VL, VR, i,j,k, DIR_X)
#else
                    ! MHD
                    call hy_uhd_HLLD(DIR_X,VL,VR,xflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM,i,j,k),speed,ierr)
                    if(ierr /= 0) call do_error("hy_uhd_HLLD", VL, VR, i,j,k, DIR_X)
#endif
                 endif
#endif
              endif

!!$#ifdef FLASH_UGLM_MHD
!!$           xflux(F06MAGX_FLUX,i,j,k) = scrchFaceXPtr(HY_N06_FACEXPTR_VAR,i,j,k)
!!$           xflux(F09GLMP_FLUX,i,j,k) = scrchFaceXPtr(HY_N09_FACEXPTR_VAR,i,j,k)
!!$#endif

!#ifdef COMPUTE_DT_FLUX
              if (hy_hydroComputeDtOption == 1) then
                 !! Call for dt calculation
                 call hy_uhd_setMinTimeStep(blockID,i,j,k,del(DIR_X),speed)
              endif
!#endif

              !! Artificial viscosity as in PPM, Colella and Woodward, 1984.
#ifdef BDRY_VAR
              if (U(BDRY_VAR,i,j,k) < 0.0 .and. U(BDRY_VAR,i-1,j,k) < 0.0) then
#endif
              if (hy_use_avisc) then
                 if (NDIM == 1) then
                    cvisc = hy_cvisc*max(    -(U(VELX_VAR,i,  j,  k  )-U(VELX_VAR,i-1,j,  k  )),0.)
                 elseif (NDIM == 2) then
                    cvisc = hy_cvisc*max(-(    U(VELX_VAR,i,  j,  k  )-U(VELX_VAR,i-1,j,  k  ) + &
                                         0.25*(U(VELY_VAR,i,  j+1,k  )-U(VELY_VAR,i,  j-1,k  ) + &
                                               U(VELY_VAR,i-1,j+1,k  )-U(VELY_VAR,i-1,j-1,k  ))*del(DIR_X)/del(DIR_Y)),&
                                         0.)
                 else
                    cvisc = hy_cvisc*max(-(     U(VELX_VAR,i,  j,  k  )-U(VELX_VAR,i-1,j,  k  ) + &
                                         0.25*( U(VELY_VAR,i,  j+1,k  )-U(VELY_VAR,i,  j-1,k  ) + &
                                                U(VELY_VAR,i-1,j+1,k  )-U(VELY_VAR,i-1,j-1,k  ))*del(DIR_X)/del(DIR_Y) + &
                                         0.25*( U(VELZ_VAR,i,  j,  k+1)-U(VELZ_VAR,i,  j,  k-1) + &
                                                U(VELZ_VAR,i-1,j,  k+1)-U(VELZ_VAR,i-1,j,  k-1))*del(DIR_X)/del(DIR_Z)), &
                                         0.)
                 endif
              
                 xflux(F01DENS_FLUX:F05ENER_FLUX,i,j,k) = &
                      xflux(F01DENS_FLUX:F05ENER_FLUX,i,j,k) &
                      +cvisc*(/U(DENS_VAR,i-1,j,k)                    -U(DENS_VAR,i,j,k)&
                              ,U(DENS_VAR,i-1,j,k)*U(VELX_VAR,i-1,j,k)-U(DENS_VAR,i,j,k)*U(VELX_VAR,i,j,k)&
                              ,U(DENS_VAR,i-1,j,k)*U(VELY_VAR,i-1,j,k)-U(DENS_VAR,i,j,k)*U(VELY_VAR,i,j,k)&
                              ,U(DENS_VAR,i-1,j,k)*U(VELZ_VAR,i-1,j,k)-U(DENS_VAR,i,j,k)*U(VELZ_VAR,i,j,k)&
                              ,U(DENS_VAR,i-1,j,k)*U(ENER_VAR,i-1,j,k)-U(DENS_VAR,i,j,k)*U(ENER_VAR,i,j,k)/)
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
                 xflux(F07MAGY_FLUX:F08MAGZ_FLUX,i,j,k) = &
                      xflux(F07MAGY_FLUX:F08MAGZ_FLUX,i,j,k) &
                      +cvisc*(/U(MAGY_VAR,i-1,j,k)                    -U(MAGY_VAR,i,j,k)&
                              ,U(MAGZ_VAR,i-1,j,k)                    -U(MAGZ_VAR,i,j,k) /)
#ifdef FLASH_UGLM_MHD
                 xflux(HY_GLMP_FLUX,i,j,k) = &
                      xflux(HY_GLMP_FLUX,i,j,k) &
                      +cvisc*(U(GLMP_VAR,i-1,j,k)                    -U(GLMP_VAR,i,j,k) )
#endif
#endif

              endif
#ifdef BDRY_VAR
              endif
#endif
              !! Flux for internal energy density
              !! Reference: "Simple Method to Track Pressure Accurately", S. Li, Astronum Proceeding, 2007
              !! Note that there is a typo in Li's paper related on the left and right states.
              if (hy_useAuxEintEqn) then
                 if (xflux(HY_DENS_FLUX,i,j,k) > 0.) then
                    xflux(HY_VOLU_FLUX,i,j,k) = xflux(HY_DENS_FLUX,i,j,k)/VL(HY_DENS)
                    xflux(HY_EINT_FLUX,i,j,k) = xflux(HY_DENS_FLUX,i,j,k)*VL(HY_EINT)
                 else
                    xflux(HY_VOLU_FLUX,i,j,k) = xflux(HY_DENS_FLUX,i,j,k)/VR(HY_DENS)
                    xflux(HY_EINT_FLUX,i,j,k) = xflux(HY_DENS_FLUX,i,j,k)*VR(HY_EINT)
                 endif

              endif

#ifdef FLASH_UHD_3T
              if (xflux(HY_DENS_FLUX,i,j,k) > 0.) then
                 xflux(HY_EELE_FLUX:HY_ERAD_FLUX,i,j,k) = xflux(HY_DENS_FLUX,i,j,k)*VL(HY_EELE:HY_ERAD)
              else
                 xflux(HY_EELE_FLUX:HY_ERAD_FLUX,i,j,k) = xflux(HY_DENS_FLUX,i,j,k)*VR(HY_EELE:HY_ERAD)
              endif

#endif

#if (NSPECIES+NMASS_SCALARS) > 0
              if (hy_fullSpecMsFluxHandling) then
                 do ispu = SPECIES_BEGIN, MASS_SCALARS_END !SPECIES_END
                    isph= ispu-NPROP_VARS
!!$                    print*,'i,j,k=',i,j,k,present(hy_spcL),present(hy_spcR)
!!$                    print*,'associated(hy_spcR)=',associated(hy_spcR)
!!$                    print*,'hy_SpcR(isph,i-1,j,k,DIR_X)=',i,j,k,hy_SpcR(isph,i-1,j,k,DIR_X)
                    if (xflux(HY_DENS_FLUX,i,j,k) < 0.) then
                       xflux(HY_END_FLUX+isph,i,j,k) = xflux(HY_DENS_FLUX,i,j,k)*hy_SpcL(isph,i,  j,k,DIR_X)
                    else
                       xflux(HY_END_FLUX+isph,i,j,k) = xflux(HY_DENS_FLUX,i,j,k)*hy_SpcR(isph,i-1,j,k,DIR_X)
                    endif
                 enddo
              endif
#endif


           endif ! end of if (hy_updateHydroFluxes) then

           if (hy_useDiffuse) then
              if (hy_useViscosity) then
                 call hy_uhd_addViscousFluxes&
                      (blockID,blkLimitsGC,i,j,k,xflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k),viscDynamic,DIR_X)
              endif

              if (hy_useConductivity .and. hy_addThermalFlux) then
                 call hy_uhd_addThermalFluxes&
                      (blockID,blkLimitsGC,i,j,k,xflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k),cond,DIR_X)
              endif

#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
              if (hy_useMagneticResistivity) then
                 call hy_uhd_addResistiveFluxes&
                      (blockID,blkLimitsGC,i,j,k,xflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k),magResist,DIR_X)
              endif
#endif
           endif

#ifdef FLASH_USM_MHD           
           if ((hy_useBiermann .or. hy_useBiermann1T) .and. (.not. hy_biermannSource)) then
              call hy_uhd_addBiermannBatteryTerms &
                   (blockID,blkLimitsGC,i,j,k,xflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k), DIR_X)
           endif 
#endif
           
        enddo
#if NDIM == 1
        !$omp end do
#endif
     enddo
#if NDIM == 2
     !$omp end do
#endif
  enddo
#if NDIM == 3
  !$omp end do
#endif


#ifdef FLASH_USM_MHD
  if (present(lastCall)) then
     if (lastCall) then
#endif
        if (hy_useAuxEintEqn .OR. hy_geometry /= CARTESIAN) then
           !! Obtain an averaged pressure for internal energy update in hy_uhd_unsplitUpdate
           
#if NDIM == 3
           !$omp do schedule(static)
#endif
           do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
#if NDIM == 2
              !$omp do schedule(static)
#endif
              do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
#if NDIM == 1
                 !$omp do schedule(static)
#endif
                 do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

                    !! For non-cartesian geometries
                    if (hy_geometry /= CARTESIAN) then
                       select case(hy_geometry) ! First, select whether y or z is phi-direction
                       case(CYLINDRICAL,POLAR)
                          alpha = 1.
                       case(SPHERICAL)
                          alpha = 2.
                       end select
                    endif !end of non-Cartesian support

                    if (hy_geometry /= CARTESIAN) then

                       presPlus = scrchFaceXPtr(HY_P05_FACEXPTR_VAR,i,j,k)
                       presMinus =scrchFaceXPtr(HY_N05_FACEXPTR_VAR,i,j,k)

#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
                       presPlus = presPlus + &
                            0.5*(scrchFaceXPtr(HY_P06_FACEXPTR_VAR,i,j,k)**2+&
                            scrchFaceXPtr(HY_P07_FACEXPTR_VAR,i,j,k)**2+&
                            scrchFaceXPtr(HY_P08_FACEXPTR_VAR,i,j,k)**2)

                       presMinus= presMinus + &
                            0.5*(scrchFaceXPtr(HY_N06_FACEXPTR_VAR,i,j,k)**2+&
                            scrchFaceXPtr(HY_N07_FACEXPTR_VAR,i,j,k)**2+&
                            scrchFaceXPtr(HY_N08_FACEXPTR_VAR,i,j,k)**2)
#endif
                       weightMinus = xLeft(i); weightPlus = xRight(i)
                       weightSum = 2.0*xCenter(i)
                       if (hy_geometry == SPHERICAL) then
                          weightMinus = xLeft(i) * xLeft(i); weightPlus=xRight(i)*xRight(i)
                          weightSum = weightMinus + weightPlus
                       end if

                       scrch_Ptr(HY_VAR2_SCRATCHCTR_VAR,i,j,k) = &
                            (weightMinus*presMinus + weightPlus*presPlus)/weightSum

                    else
                       scrch_Ptr(HY_VAR2_SCRATCHCTR_VAR,i,j,k) = &
                            0.5*( scrchFaceXPtr(HY_N05_FACEXPTR_VAR,i,j,k) &
                                 +scrchFaceXPtr(HY_P05_FACEXPTR_VAR,i,j,k))
                    endif

                 enddo
#if NDIM == 1
                 !$omp end do
#endif
              enddo
#if NDIM == 2
              !$omp end do
#endif
           enddo
#if NDIM == 3
           !$omp end do
#endif
        endif ! end of if (hy_useAuxEintEqn)

#ifndef FLASH_GRID_UG
        if (hy_fullRiemannStateArrays) then

        !use scrch arrays for flux conservation
           !$omp workshare
           scrchFaceXPtr(HY_P01_FACEXPTR_VAR:HY_P01_FACEXPTR_VAR+HY_SCRATCH_NUM-1,&
                                             i0:imax+1,j0:jmax,k0:kmax) &
            = xflux(HY_DENS_FLUX:HY_END_FLUX,i0:imax+1,j0:jmax,k0:kmax)
           !$omp end workshare
        end if ! (hy_fullRiemannStateArrays)
#endif


#ifdef FLASH_USM_MHD
     endif !end of if (lastCall)
  endif !end of if-present(lastCall)
#endif



#if NDIM >= 2

  !! Calculate y-flux
  !$omp single
  normDir=0; normDir(DIR_Y)=1
  tranDir=2; tranDir(DIR_Y)=0
  tranDir=tranDir*kUSM
  tranDir(2)=tranDir(2)*k2
  tranDir(3)=tranDir(3)*k3
  !$omp end single 

#if NDIM == 3
  !$omp do schedule(static)
#endif
  do k=blkLimits(LOW,KAXIS)-tranDir(DIR_Z),blkLimits(HIGH,KAXIS)+tranDir(DIR_Z)+normDir(DIR_Z)
#if NDIM == 2
     !$omp do schedule(static)
#endif
     do j=blkLimits(LOW,JAXIS)-tranDir(DIR_Y),blkLimits(HIGH,JAXIS)+tranDir(DIR_Y)+normDir(DIR_Y)
        do i=blkLimits(LOW,IAXIS)-tranDir(DIR_X),blkLimits(HIGH,IAXIS)+tranDir(DIR_X)+normDir(DIR_X)

        if (hy_updateHydroFluxes) then

           VL(HY_DENS:HY_END_VARS-kGrav)=scrchFaceYPtr(HY_P01_FACEYPTR_VAR:HY_P01_FACEYPTR_VAR+HY_SCRATCH_NUM-1,i,j-1,k)
           VR(HY_DENS:HY_END_VARS-kGrav)=scrchFaceYPtr(HY_N01_FACEYPTR_VAR:HY_N01_FACEYPTR_VAR+HY_SCRATCH_NUM-1,i,j,  k)

#ifdef BDRY_VAR
           ! solid internal boundary
           if (U(BDRY_VAR,i,j,k) > 0.0 .and. U(BDRY_VAR,i,j-1,k) < 0.0) then
              VR(HY_DENS:HY_END_VARS-kGrav) = VL(HY_DENS:HY_END_VARS-kGrav)
              VR(HY_VELY) = -VL(HY_VELY)
           end if

           if (U(BDRY_VAR,i,j,k) < 0.0 .and. U(BDRY_VAR,i,j-1,k) > 0.0) then
              VL(HY_DENS:HY_END_VARS-kGrav) = VR(HY_DENS:HY_END_VARS-kGrav)
              VL(HY_VELY) = -VR(HY_VELY)
           end if
#endif

           if (hy_EOSforRiemann) then
              ! Construct initial guess for temperature
              t_start = 0.5*(U(TEMP_VAR,i,j-1,k) + U(TEMP_VAR,i,j,k))
              if (hy_fullSpecMsFluxHandling .AND. hy_numXN > 0) then
                 call faceStatesEos(VL,t_start,hy_SpcL(:,i,j  ,k,DIR_Y))
                 call faceStatesEos(VR,t_start,hy_SpcR(:,i,j-1,k,DIR_Y))
              else
                 call faceStatesEos(VL,t_start)
                 call faceStatesEos(VR,t_start)
              end if
           end if

           if (hy_RiemannSolver == ROE) then
              call hy_uhd_Roe(DIR_Y,VL,VR,yflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM,i,j,k),speed,ierr)
              if(ierr /= 0) call do_error("hy_uhd_Roe", VL, VR, i,j,k, DIR_Y)

           elseif (hy_RiemannSolver == HLL) then
              call hy_uhd_HLL(DIR_Y,VL,VR,yflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM,i,j,k),speed,ierr)
              if(ierr /= 0) call do_error("hy_uhd_HLL", VL, VR, i,j,k, DIR_Y)

           elseif (hy_RiemannSolver == HLLC) then
              call hy_uhd_HLLC(DIR_Y,VL,VR,yflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM,i,j,k),speed,ierr)
              if(ierr /= 0) call do_error("hy_uhd_HLLC", VL, VR, i,j,k, DIR_Y)

#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
           elseif (hy_RiemannSolver == HLLD) then
              call hy_uhd_HLLD(DIR_Y,VL,VR,yflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM,i,j,k),speed,ierr)
              if(ierr /= 0) call do_error("hy_uhd_HLLD", VL, VR, i,j,k, DIR_Y)
#endif
           elseif (hy_RiemannSolver == MARQ) then
              call hy_uhd_Marquina(DIR_Y,VL,VR,yflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM,i,j,k),speed,ierr)
              if(ierr /= 0) call do_error("hy_uhd_Marquina", VL, VR, i,j,k, DIR_Y)

           elseif (hy_RiemannSolver == MARM) then
              call hy_uhd_MarquinaModified(DIR_Y,VL,VR,yflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM,i,j,k),speed,ierr)
              if(ierr /= 0) call do_error("hy_uhd_Marquina", VL, VR, i,j,k, DIR_Y)

           elseif (hy_RiemannSolver == LLF) then
              call hy_uhd_LLF(DIR_Y,VL,VR,yflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM,i,j,k),speed,ierr)
              if(ierr /= 0) call do_error("hy_uhd_LLF", VL, VR, i,j,k, DIR_Y)

#ifdef SHOK_VAR
           elseif (hy_RiemannSolver == HYBR) then
              if (U(SHOK_VAR,i,j-1,k) + U(SHOK_VAR,i,j,k) > 0.) then
                 !! use diffusive HLL solver for a local strong shock/rarefaction region
                 call hy_uhd_HLL(DIR_Y,VL,VR,yflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM,i,j,k),speed,ierr)
                 if(ierr /= 0) call do_error("hy_uhd_HLL", VL, VR, i,j,k, DIR_Y)
              else
#ifdef FLASH_UHD_HYDRO
                 ! for hydro
                 call hy_uhd_HLLC(DIR_Y,VL,VR,yflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM,i,j,k),speed,ierr)
                 if(ierr /= 0) call do_error("hy_uhd_HLLC", VL, VR, i,j,k, DIR_Y)
#else
                 ! for MHD
                 call hy_uhd_HLLD(DIR_Y,VL,VR,yflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM,i,j,k),speed,ierr)
                 if(ierr /= 0) call do_error("hy_uhd_HLLD", VL, VR, i,j,k, DIR_Y)
#endif
              endif
#endif
           endif


!!$#ifdef FLASH_UGLM_MHD
!!$           yflux(F07MAGY_FLUX,i,j,k) = scrchFaceYPtr(HY_N07_FACEYPTR_VAR,i,j,k)
!!$           yflux(F09GLMP_FLUX,i,j,k) = scrchFaceYPtr(HY_N09_FACEYPTR_VAR,i,j,k)
!!$#endif

!#ifdef COMPUTE_DT_FLUX
           if (hy_hydroComputeDtOption == 1) then
              !! Call for dt calculation
              if (hy_geometry == CARTESIAN .OR. hy_geometry == CYLINDRICAL) then 
                 ! the 'y' coordinate is not angular in Cartesian and cylindrical coords
                 dy = del(DIR_Y)
              else ! Angular coordinates in 2D: Spherical or Polar
                 dy = xCenter(i)*del(DIR_Y)
              endif
              call hy_uhd_setMinTimeStep(blockID,i,j,k,dy,speed)
           endif
!#endif

           !! Artificial viscosity as in PPM, Colella and Woodward, 1984.
#ifdef BDRY_VAR
           if (U(BDRY_VAR,i,j,k) < 0.0 .and. U(BDRY_VAR,i,j-1,k) < 0.0) then
#endif
           if (hy_use_avisc) then
              if (NDIM == 2) then
                 cvisc = hy_cvisc*max(-(0.25*( U(VELX_VAR,i+1,j,  k  )-U(VELX_VAR,i-1,j,  k  ) + &
                                               U(VELX_VAR,i+1,j-1,k  )-U(VELX_VAR,i-1,j-1,k  ))*del(DIR_Y)/del(DIR_X) + &
                                               U(VELY_VAR,i,  j,  k  )-U(VELY_VAR,i,  j-1,k  )), &
                                      0.)
              else
                 cvisc = hy_cvisc*max(-(0.25*( U(VELX_VAR,i+1,j,  k  )-U(VELX_VAR,i-1,j,  k  ) + &
                                               U(VELX_VAR,i+1,j-1,k  )-U(VELX_VAR,i-1,j-1,k  ))*del(DIR_Y)/del(DIR_X) + &
                                               U(VELY_VAR,i,  j,  k  )-U(VELY_VAR,i,  j-1,k  ) + &
                                        0.25*( U(VELZ_VAR,i,  j,  k+1)-U(VELZ_VAR,i,  j,  k-1) + &
                                               U(VELZ_VAR,i,  j-1,k+1)-U(VELZ_VAR,i,  j-1,k-1))*del(DIR_Y)/del(DIR_Z)), &
                                      0.)

              endif

              yflux(F01DENS_FLUX:F05ENER_FLUX,i,j,k) = &
                   yflux(F01DENS_FLUX:F05ENER_FLUX,i,j,k) &
                   +cvisc*(/U(DENS_VAR,i,j-1,k)                    -U(DENS_VAR,i,j,k)&
                           ,U(DENS_VAR,i,j-1,k)*U(VELX_VAR,i,j-1,k)-U(DENS_VAR,i,j,k)*U(VELX_VAR,i,j,k)&
                           ,U(DENS_VAR,i,j-1,k)*U(VELY_VAR,i,j-1,k)-U(DENS_VAR,i,j,k)*U(VELY_VAR,i,j,k)&
                           ,U(DENS_VAR,i,j-1,k)*U(VELZ_VAR,i,j-1,k)-U(DENS_VAR,i,j,k)*U(VELZ_VAR,i,j,k)&
                           ,U(DENS_VAR,i,j-1,k)*U(ENER_VAR,i,j-1,k)-U(DENS_VAR,i,j,k)*U(ENER_VAR,i,j,k)/)
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
              yflux(F06MAGX_FLUX:F08MAGZ_FLUX,i,j,k) = &
                   yflux(F06MAGX_FLUX:F08MAGZ_FLUX,i,j,k) &
                   +cvisc*(/U(MAGX_VAR,i,j-1,k)                    -U(MAGX_VAR,i,j,k)&
                           ,0.&
                           ,U(MAGZ_VAR,i,j-1,k)                    -U(MAGZ_VAR,i,j,k) /)
#ifdef FLASH_UGLM_MHD
              yflux(HY_GLMP_FLUX,i,j,k) = &
                   yflux(HY_GLMP_FLUX,i,j,k) &
                   +cvisc*(U(GLMP_VAR,i,j-1,k)                    -U(GLMP_VAR,i,j,k) )
#endif
#endif

           endif
#ifdef BDRY_VAR
           endif
#endif

           if (hy_useAuxEintEqn) then
              !! Flux for internal energy density
              !! Reference: "Simple Method to Track Pressure Accurately", S. Li, Astronum Proceeding, 2007
              if (yflux(HY_DENS_FLUX,i,j,k) > 0.) then
                 yflux(HY_VOLU_FLUX,i,j,k) = yflux(HY_DENS_FLUX,i,j,k)/VL(HY_DENS)
                 yflux(HY_EINT_FLUX,i,j,k) = yflux(HY_DENS_FLUX,i,j,k)*VL(HY_EINT)
              else
                 yflux(HY_VOLU_FLUX,i,j,k) = yflux(HY_DENS_FLUX,i,j,k)/VR(HY_DENS)
                 yflux(HY_EINT_FLUX,i,j,k) = yflux(HY_DENS_FLUX,i,j,k)*VR(HY_EINT)
              endif
           endif

#ifdef FLASH_UHD_3T
           if (yflux(HY_DENS_FLUX,i,j,k) > 0.) then
              yflux(HY_EELE_FLUX:HY_ERAD_FLUX,i,j,k) = yflux(HY_DENS_FLUX,i,j,k)*VL(HY_EELE:HY_ERAD)
           else
              yflux(HY_EELE_FLUX:HY_ERAD_FLUX,i,j,k) = yflux(HY_DENS_FLUX,i,j,k)*VR(HY_EELE:HY_ERAD)
           endif
#endif

#if (NSPECIES+NMASS_SCALARS) > 0
           if (hy_fullSpecMsFluxHandling) then
              do ispu = SPECIES_BEGIN, MASS_SCALARS_END !SPECIES_END
                 isph= ispu-NPROP_VARS
                 if (yflux(HY_DENS_FLUX,i,j,k) < 0.) then
                    yflux(HY_END_FLUX+isph,i,j,k) = yflux(HY_DENS_FLUX,i,j,k)*hy_SpcL(isph,i,j,  k,DIR_Y)
                 else
                    yflux(HY_END_FLUX+isph,i,j,k) = yflux(HY_DENS_FLUX,i,j,k)*hy_SpcR(isph,i,j-1,k,DIR_Y)
                 endif
              enddo
           endif
#endif

        endif !end of if (hy_updateHydroFluxes) then

           if (hy_useDiffuse) then
              if (hy_useViscosity) then
                 call hy_uhd_addViscousFluxes&
                      (blockID,blkLimitsGC,i,j,k,yflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k),viscDynamic,DIR_Y)
              endif

              if (hy_useConductivity .and. hy_addThermalFlux) then
                 call hy_uhd_addThermalFluxes&
                      (blockID,blkLimitsGC,i,j,k,yflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k),cond,DIR_Y)
              endif

#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
              if (hy_useMagneticResistivity) then
                 call hy_uhd_addResistiveFluxes&
                      (blockID,blkLimitsGC,i,j,k,yflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k),magResist,DIR_Y)
              endif
#endif
           endif

#ifdef FLASH_USM_MHD           
           if ((hy_useBiermann .or. hy_useBiermann1T) .and. (.not. hy_biermannSource)) then
              call hy_uhd_addBiermannBatteryTerms &
                   (blockID,blkLimitsGC,i,j,k,yflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k), DIR_Y)
           endif 
#endif
           
        enddo
     enddo
#if NDIM == 2
     !$omp end do
#endif
  enddo
#if NDIM == 3
  !$omp end do
#endif



#ifdef FLASH_USM_MHD
  if (present(lastCall)) then
     if (lastCall) then
#endif
        if (hy_useAuxEintEqn) then
           !! Average out pressures at j+1/2 (P05_FACEY) and j-1/2 (N05_FACEY) and store them to
           !! scrchFaceYPtr(HY_N05_FACEYPTR_VAR,i,j,k) before the P*_FACEY array is potentially wiped out and
           !! potentially reused for storing fluxes.
           !$omp workshare
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
           !! First, increase the pressures at j+1/2 and j-1/2 by Y-face-centered magnetic pressures
           scrchFaceYPtr(HY_P05_FACEYPTR_VAR,i0:imax,j0:jmax,k0:kmax) = &
                    scrchFaceYPtr(HY_P05_FACEYPTR_VAR,i0:imax,j0:jmax,k0:kmax)  &
                    + 0.5*(scrchFaceYPtr(HY_P06_FACEYPTR_VAR,i0:imax,j0:jmax,k0:kmax)**2+&
                           scrchFaceYPtr(HY_P07_FACEYPTR_VAR,i0:imax,j0:jmax,k0:kmax)**2+&
                           scrchFaceYPtr(HY_P08_FACEYPTR_VAR,i0:imax,j0:jmax,k0:kmax)**2)
           scrchFaceYPtr(HY_N05_FACEYPTR_VAR,i0:imax,j0:jmax,k0:kmax) = &
                    scrchFaceYPtr(HY_N05_FACEYPTR_VAR,i0:imax,j0:jmax,k0:kmax)  &
                    + 0.5*(scrchFaceYPtr(HY_N06_FACEYPTR_VAR,i0:imax,j0:jmax,k0:kmax)**2+&
                           scrchFaceYPtr(HY_N07_FACEYPTR_VAR,i0:imax,j0:jmax,k0:kmax)**2+&
                           scrchFaceYPtr(HY_N08_FACEYPTR_VAR,i0:imax,j0:jmax,k0:kmax)**2)
#endif
           scrchFaceYPtr(HY_N05_FACEYPTR_VAR,:,:,:) = &
                (scrchFaceYPtr(HY_N05_FACEYPTR_VAR,:,:,:) +  scrchFaceYPtr(HY_P05_FACEYPTR_VAR,:,:,:))*0.5
           !$omp end workshare
        endif !end of if (hy_useAuxEintEqn) then
#ifndef FLASH_GRID_UG
        if (hy_fullRiemannStateArrays) then
!keeping the option to use scrch arrays for flux conservation
           !$omp workshare
        !scrchFaceYPtr(HY_P01_FACEYPTR_VAR:HY_P01_FACEYPTR_VAR+HY_SCRATCH_NUM-1,:,:,:) = 0.
           !! re-use YP array for storing fluxes
           scrchFaceYPtr(HY_P01_FACEYPTR_VAR:HY_P01_FACEYPTR_VAR+HY_SCRATCH_NUM-1,&
                                             i0:imax,jbeg:jend+1,k0:kmax)   &
            = yflux(HY_DENS_FLUX:HY_END_FLUX,i0:imax,jbeg:jend+1,k0:kmax)
           !$omp end workshare
        end if ! (hy_fullRiemannStateArrays)
#endif                

#ifdef FLASH_USM_MHD
     endif !end of if (lastCall)
  endif !end of if-present(lastCall)
#endif


#if NDIM == 3
  !! Calculate z-flux
  !$omp single
  normDir=0; normDir(DIR_Z)=1
  tranDir=2; tranDir(DIR_Z)=0
  tranDir=tranDir*kUSM
  tranDir(2)=tranDir(2)*k2
  tranDir(3)=tranDir(3)*k3
  !$omp end single
  
  !$omp do schedule(static)
  do k=blkLimits(LOW,KAXIS)-tranDir(DIR_Z),blkLimits(HIGH,KAXIS)+tranDir(DIR_Z)+normDir(DIR_Z)
     do j=blkLimits(LOW,JAXIS)-tranDir(DIR_Y),blkLimits(HIGH,JAXIS)+tranDir(DIR_Y)+normDir(DIR_Y)
        do i=blkLimits(LOW,IAXIS)-tranDir(DIR_X),blkLimits(HIGH,IAXIS)+tranDir(DIR_X)+normDir(DIR_X)

        if (hy_updateHydroFluxes) then

           VL(HY_DENS:HY_END_VARS-kGrav)=scrchFaceZPtr(HY_P01_FACEZPTR_VAR:HY_P01_FACEZPTR_VAR+HY_SCRATCH_NUM-1,i,j,k-1)
           VR(HY_DENS:HY_END_VARS-kGrav)=scrchFaceZPtr(HY_N01_FACEZPTR_VAR:HY_N01_FACEZPTR_VAR+HY_SCRATCH_NUM-1,i,j,k)

#ifdef BDRY_VAR
           ! solid internal boundary
           if (U(BDRY_VAR,i,j,k) > 0.0 .and. U(BDRY_VAR,i,j,k-1) < 0.0) then
              VR(HY_DENS:HY_END_VARS-kGrav) = VL(HY_DENS:HY_END_VARS-kGrav)
              VR(HY_VELZ) = -VL(HY_VELZ)
           end if

           if (U(BDRY_VAR,i,j,k) < 0.0 .and. U(BDRY_VAR,i,j,k-1) > 0.0) then
              VL(HY_DENS:HY_END_VARS-kGrav) = VR(HY_DENS:HY_END_VARS-kGrav)
              VL(HY_VELZ) = -VR(HY_VELZ)
           end if
#endif
           if (hy_EOSforRiemann) then
              ! Construct initial guess for temperature
              t_start = 0.5*(U(TEMP_VAR,i,j,k-1) + U(TEMP_VAR,i,j,k))
              if (hy_fullSpecMsFluxHandling .AND. hy_numXN > 0) then
                 call faceStatesEos(VL,t_start,hy_SpcL(:,i,j,k  ,DIR_Z))
                 call faceStatesEos(VR,t_start,hy_SpcR(:,i,j,k-1,DIR_Z))
              else
                 call faceStatesEos(VL,t_start)
                 call faceStatesEos(VR,t_start)
              end if
           end if

           if (hy_RiemannSolver == ROE) then
              call hy_uhd_Roe(DIR_Z,VL,VR,zflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM,i,j,k),speed,ierr)
              if(ierr /= 0) call do_error("hy_uhd_Roe", VL, VR, i,j,k, DIR_Z)

           elseif (hy_RiemannSolver == HLL) then
              call hy_uhd_HLL(DIR_Z,VL,VR,zflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM,i,j,k),speed,ierr)
              if(ierr /= 0) call do_error("hy_uhd_HLL", VL, VR, i,j,k, DIR_Z)

           elseif (hy_RiemannSolver == HLLC) then
              call hy_uhd_HLLC(DIR_Z,VL,VR,zflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM,i,j,k),speed,ierr)
              if(ierr /= 0) call do_error("hy_uhd_HLLC", VL, VR, i,j,k, DIR_Z)

#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
           elseif (hy_RiemannSolver == HLLD) then
              call hy_uhd_HLLD(DIR_Z,VL,VR,zflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM,i,j,k),speed,ierr)
              if(ierr /= 0) call do_error("hy_uhd_HLLD", VL, VR, i,j,k, DIR_Z)
#endif
           elseif (hy_RiemannSolver == MARQ) then
              call hy_uhd_Marquina(DIR_Z,VL,VR,zflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM,i,j,k),speed,ierr)
              if(ierr /= 0) call do_error("hy_uhd_Marquina", VL, VR, i,j,k, DIR_Z)

           elseif (hy_RiemannSolver == MARM) then
              call hy_uhd_MarquinaModified(DIR_Z,VL,VR,zflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM,i,j,k),speed,ierr)
              if(ierr /= 0) call do_error("hy_uhd_Marquina", VL, VR, i,j,k, DIR_Z)

           elseif (hy_RiemannSolver == LLF) then
              call hy_uhd_LLF(DIR_Z,VL,VR,zflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM,i,j,k),speed,ierr)
              if(ierr /= 0) call do_error("hy_uhd_LLF", VL, VR, i,j,k, DIR_Z)

#ifdef SHOK_VAR
           elseif (hy_RiemannSolver == HYBR) then
              if (U(SHOK_VAR,i,j,k-1) + U(SHOK_VAR,i,j,k) > 0.) then
                 !! use diffusive HLL solver for a local strong shock/rarefaction region
                 call hy_uhd_HLL(DIR_Z,VL,VR,zflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM,i,j,k),speed,ierr)
                 if(ierr /= 0) call do_error("hy_uhd_HLL", VL, VR, i,j,k, DIR_Z)
              else
#ifdef FLASH_UHD_HYDRO
                 ! for hydro
                 call hy_uhd_HLLC(DIR_Z,VL,VR,zflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM,i,j,k),speed,ierr)
                 if(ierr /= 0) call do_error("hy_uhd_HLLC", VL, VR, i,j,k, DIR_Z)
#else
                 ! for MHD
                 call hy_uhd_HLLD(DIR_Z,VL,VR,zflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM,i,j,k),speed,ierr)
                 if(ierr /= 0) call do_error("hy_uhd_HLLD", VL, VR, i,j,k, DIR_Z)
#endif
              endif
#endif

           endif
 
!#ifdef COMPUTE_DT_FLUX
           if (hy_hydroComputeDtOption == 1) then
              !! Call for dt calculation
              if (hy_geometry == CARTESIAN) then 
                 ! the 'y' coordinate is not angular in Cartesian and cylindrical coords
                 dz = del(DIR_Z)
              elseif (hy_geometry == CYLINDRICAL) then
                 dz = xCenter(i)*del(DIR_Z)
              else ! Angular coordinates in 2D: Spherical or Polar
                 dz = xCenter(i)*sin(yCenter(j))*del(DIR_Z) ! z is phi
              endif
              call hy_uhd_setMinTimeStep(blockID,i,j,k,dz,speed)
           endif
!#endif
           !! Artificial viscosity as in PPM, Colella and Woodward, 1984.
#ifdef BDRY_VAR
           if (U(BDRY_VAR,i,j,k) < 0.0 .and. U(BDRY_VAR,i,j,k-1) < 0.0) then
#endif
           if (hy_use_avisc) then
              cvisc = hy_cvisc*max(-(0.25*( U(VELX_VAR,i+1,j,  k  )-U(VELX_VAR,i-1,j,  k  ) + &
                                            U(VELX_VAR,i+1,j,  k-1)-U(VELX_VAR,i-1,j,  k-1))*del(DIR_Z)/del(DIR_X) + &
                                     0.25*( U(VELY_VAR,i,  j+1,k  )-U(VELY_VAR,i,  j-1,k  ) + &
                                            U(VELY_VAR,i,  j+1,k-1)-U(VELY_VAR,i,  j-1,k-1))*del(DIR_Z)/del(DIR_Y) + &
                                            U(VELZ_VAR,i,  j,  k  )-U(VELZ_VAR,i,  j,  k-1)),&
                                   0.)

              zflux(F01DENS_FLUX:F05ENER_FLUX,i,j,k) = &
                   zflux(F01DENS_FLUX:F05ENER_FLUX,i,j,k) &
                   +cvisc*(/U(DENS_VAR,i,j,k-1)                    -U(DENS_VAR,i,j,k)&
                           ,U(DENS_VAR,i,j,k-1)*U(VELX_VAR,i,j,k-1)-U(DENS_VAR,i,j,k)*U(VELX_VAR,i,j,k)&
                           ,U(DENS_VAR,i,j,k-1)*U(VELY_VAR,i,j,k-1)-U(DENS_VAR,i,j,k)*U(VELY_VAR,i,j,k)&
                           ,U(DENS_VAR,i,j,k-1)*U(VELZ_VAR,i,j,k-1)-U(DENS_VAR,i,j,k)*U(VELZ_VAR,i,j,k)&
                           ,U(DENS_VAR,i,j,k-1)*U(ENER_VAR,i,j,k-1)-U(DENS_VAR,i,j,k)*U(ENER_VAR,i,j,k)/)
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
              zflux(F06MAGX_FLUX:F07MAGY_FLUX,i,j,k) = &
                   zflux(F06MAGX_FLUX:F07MAGY_FLUX,i,j,k) &
                   +cvisc*(/U(MAGX_VAR,i,j,k-1)                    -U(MAGX_VAR,i,j,k)&
                           ,U(MAGY_VAR,i,j,k-1)                    -U(MAGY_VAR,i,j,k) /)
#ifdef FLASH_UGLM_MHD
              zflux(HY_GLMP_FLUX,i,j,k) = &
                   zflux(HY_GLMP_FLUX,i,j,k) &
                   +cvisc*(U(GLMP_VAR,i,j,k-1)                    -U(GLMP_VAR,i,j,k) )
#endif
#endif
           endif
#ifdef BDRY_VAR
           endif
#endif

           if (hy_useAuxEintEqn) then
              !! Flux for internal energy density
              !! Reference: "Simple Method to Track Pressure Accurately", S. Li, Astronum Proceeding, 2007
              if (zflux(HY_DENS_FLUX,i,j,k) > 0.) then
                 zflux(HY_VOLU_FLUX,i,j,k) = zflux(HY_DENS_FLUX,i,j,k)/VL(HY_DENS)
                 zflux(HY_EINT_FLUX,i,j,k) = zflux(HY_DENS_FLUX,i,j,k)*VL(HY_EINT)
              else
                 zflux(HY_VOLU_FLUX,i,j,k) = zflux(HY_DENS_FLUX,i,j,k)/VR(HY_DENS)
                 zflux(HY_EINT_FLUX,i,j,k) = zflux(HY_DENS_FLUX,i,j,k)*VR(HY_EINT)
              endif
           endif

#ifdef FLASH_UHD_3T
           if (zflux(HY_DENS_FLUX,i,j,k) > 0.) then
              zflux(HY_EELE_FLUX:HY_ERAD_FLUX,i,j,k) = zflux(HY_DENS_FLUX,i,j,k)*VL(HY_EELE:HY_ERAD)
           else
              zflux(HY_EELE_FLUX:HY_ERAD_FLUX,i,j,k) = zflux(HY_DENS_FLUX,i,j,k)*VR(HY_EELE:HY_ERAD)
           endif
#endif


#if (NSPECIES+NMASS_SCALARS) > 0
           if (hy_fullSpecMsFluxHandling) then
              do ispu = SPECIES_BEGIN, MASS_SCALARS_END !SPECIES_END
                 isph= ispu-NPROP_VARS
                 if (zflux(HY_DENS_FLUX,i,j,k) < 0.) then
                    zflux(HY_END_FLUX+isph,i,j,k) = zflux(HY_DENS_FLUX,i,j,k)*hy_SpcL(isph,i,j,k,  DIR_Z)
                 else
                    zflux(HY_END_FLUX+isph,i,j,k) = zflux(HY_DENS_FLUX,i,j,k)*hy_SpcR(isph,i,j,k-1,DIR_Z)
                 endif
              enddo
           endif
#endif

        endif !end of if (hy_updateHydroFluxes) then

           if (hy_useDiffuse) then
              if (hy_useViscosity) then
                 call hy_uhd_addViscousFluxes&
                      (blockID,blkLimitsGC,i,j,k,zflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k),viscDynamic,DIR_Z)
              endif

              if (hy_useConductivity .and. hy_addThermalFlux) then
                 call hy_uhd_addThermalFluxes&
                      (blockID,blkLimitsGC,i,j,k,zflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k),cond,DIR_Z)
              endif

#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
              if (hy_useMagneticResistivity) then
                 call hy_uhd_addResistiveFluxes&
                      (blockID,blkLimitsGC,i,j,k,zflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k),magResist,DIR_Z)
              endif
#endif
           endif

#ifdef FLASH_USM_MHD
           
           if ((hy_useBiermann .or. hy_useBiermann1T) .and. (.not. hy_biermannSource)) then
              call hy_uhd_addBiermannBatteryTerms &
                   (blockID,blkLimitsGC,i,j,k,zflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k), DIR_Z)
           endif 
#endif

        enddo
     enddo
  enddo
  !$omp end do


#ifdef FLASH_USM_MHD
  if (present(lastCall)) then
     if (lastCall) then
#endif
        if (hy_useAuxEintEqn) then
           !! Average out pressures at k+1/2 (P05_FACEZ) and k-1/2 (N05_FACEZ) and store them to
           !! scrchFaceZPtr(HY_N05_FACEZPTR_VAR,i,j,k) before the P*_FACEZ array is potentially wiped out and
           !! reused for storing fluxes.
           !$omp workshare
#ifdef FLASH_USM_MHD
           scrchFaceZPtr(HY_P05_FACEZPTR_VAR,i0:imax,j0:jmax,k0:kmax) = &
                    scrchFaceZPtr(HY_P05_FACEZPTR_VAR,i0:imax,j0:jmax,k0:kmax)  &
                    + 0.5*(scrchFaceZPtr(HY_P06_FACEZPTR_VAR,i0:imax,j0:jmax,k0:kmax)**2+&
                           scrchFaceZPtr(HY_P07_FACEZPTR_VAR,i0:imax,j0:jmax,k0:kmax)**2+&
                           scrchFaceZPtr(HY_P08_FACEZPTR_VAR,i0:imax,j0:jmax,k0:kmax)**2)
           scrchFaceZPtr(HY_N05_FACEZPTR_VAR,i0:imax,j0:jmax,k0:kmax) = &
                    scrchFaceZPtr(HY_N05_FACEZPTR_VAR,i0:imax,j0:jmax,k0:kmax)  &
                    + 0.5*(scrchFaceZPtr(HY_N06_FACEZPTR_VAR,i0:imax,j0:jmax,k0:kmax)**2+&
                           scrchFaceZPtr(HY_N07_FACEZPTR_VAR,i0:imax,j0:jmax,k0:kmax)**2+&
                           scrchFaceZPtr(HY_N08_FACEZPTR_VAR,i0:imax,j0:jmax,k0:kmax)**2)
#endif
           scrchFaceZPtr(HY_N05_FACEZPTR_VAR,:,:,:) = &
                (scrchFaceZPtr(HY_N05_FACEZPTR_VAR,:,:,:) +  scrchFaceZPtr(HY_P05_FACEZPTR_VAR,:,:,:))*0.5
           !$omp end workshare
        endif !end of if (hy_useAuxEintEqn)

#ifndef FLASH_GRID_UG
        if (hy_fullRiemannStateArrays) then
!keeping the option to use scrch arrays for flux conservation
           !$omp workshare
           !! initialize with zero
           !scrchFaceZPtr(HY_P01_FACEZPTR_VAR:HY_P01_FACEZPTR_VAR+HY_SCRATCH_NUM-1,:,:,:) = 0.
           !! re-use ZP array for storing fluxes
           scrchFaceZPtr(HY_P01_FACEZPTR_VAR:HY_P01_FACEZPTR_VAR+HY_SCRATCH_NUM-1,&
                                             i0:imax,j0:jmax,kbeg:kend+1)   &
            = zflux(HY_DENS_FLUX:HY_END_FLUX,i0:imax,j0:jmax,kbeg:kend+1)
           !$omp end workshare
        end if ! (hy_fullRiemannStateArrays)
#endif

      
#ifdef FLASH_USM_MHD
     endif !end of if (lastCall)
  endif !end of if-present(lastCall)
#endif


#endif /* end of #if NDIM==3 */
#endif /* end of #if NDIM >= 2 */


  if (hy_useAuxEintEqn) then
  !! Average out the n+1/2 pressures at the cell interfaces to the cell center on each cell
  !! and store them to scrch_Ptr(HY_VAR1_SCRATCHCTR_VAR,:,:,:). 
  !! These averaged pressures at n+1/2 are used later for the internal energy update.
  !! Note that scrch_Ptr(HY_VAR2_SCRATCHCTR_VAR,:,:,:) holds the averaged pressures 
  !! in radial R-direction only for non-Cartesian cases for geometric source term 
  !! and should not be modified here.

  !! (1) store the averaged pressure in 'R-direction' to VAR2, and
  !! (2) store the averaged pressure in 'ALL directions' to VAR1
  !$omp workshare
  scrch_Ptr(HY_VAR1_SCRATCHCTR_VAR,:,:,:) = scrch_Ptr(HY_VAR2_SCRATCHCTR_VAR,:,:,:)
#if NDIM > 1
  scrch_Ptr(HY_VAR1_SCRATCHCTR_VAR,&
             blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) = &
       (scrch_Ptr(HY_VAR1_SCRATCHCTR_VAR,&
             blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))&
        +scrchFaceYPtr(HY_N05_FACEYPTR_VAR,&
             blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)))*0.5
#if NDIM > 2
  scrch_Ptr(HY_VAR1_SCRATCHCTR_VAR,&
             blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) = &
      (2.*scrch_Ptr(HY_VAR1_SCRATCHCTR_VAR,&
             blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))&
       +scrchFaceZPtr(HY_N05_FACEZPTR_VAR,&
             blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)))/3.
#endif
#endif
  !$omp end workshare
  endif ! end of if (hy_useAuxEintEqn) then
  !$omp end parallel
  
  !! Release pointer
  call Grid_releaseBlkPtr(blockID,U,CENTER)


contains

  subroutine do_error(solver, VL, VR, i,j,k, dir)
    use Grid_interface, ONLY: Grid_getCellCoords
    use Driver_interface, ONLY: Driver_abortFlash
    implicit none
    
    character(len=*), intent(in) :: solver
    real, intent(IN) :: VL(:)
    real, intent(IN) :: VR(:)
    integer, intent(IN) :: i, j, k, dir

    print *, "LEFT STATE:"
    print *, "DENS: ", VL(HY_DENS)
    print *, "PRES: ", VL(HY_PRES)
    print *, "GAMC: ", VL(HY_GAMC)
    print *, "VELX: ", VL(HY_VELX)
    print *, "VELY: ", VL(HY_VELY)
    print *, "VELZ: ", VL(HY_VELZ)

    print *, ""
    print *, "RIGHT STATE:"
    print *, "DENS: ", VR(HY_DENS)
    print *, "PRES: ", VR(HY_PRES)
    print *, "GAMC: ", VR(HY_GAMC)
    print *, "VELX: ", VR(HY_VELX)
    print *, "VELY: ", VR(HY_VELY)
    print *, "VELZ: ", VR(HY_VELZ)

    print *, ""
    print *, "INTERFACE: ", i, j, k

    ! Note: We only allocate these arrays here, which is fine because
    !       the code is to be aborted immediately once this routine is called.
    allocate(xcent(blkLimitsGC(HIGH, IAXIS)))
    allocate(ycent(blkLimitsGC(HIGH, JAXIS)))
    allocate(zcent(blkLimitsGC(HIGH, KAXIS)))

    call Grid_getCellCoords(IAXIS, blockId, CENTER, .true., xcent, blkLimitsGC(HIGH, IAXIS)) 
    call Grid_getCellCoords(JAXIS, blockId, CENTER, .true., ycent, blkLimitsGC(HIGH, JAXIS))
    call Grid_getCellCoords(KAXIS, blockId, CENTER, .true., zcent, blkLimitsGC(HIGH, KAXIS))

    print *, "NEIGBORING CELLS:"

    if(dir == DIR_X) then
       print *, "X DIRECTION"
       call write_cell_info(i-3,j,k)
       call write_cell_info(i-2,j,k)
       call write_cell_info(i-1,j,k)
       call write_cell_info(i+0,j,k)
       call write_cell_info(i+1,j,k)
       call write_cell_info(i+2,j,k)
    elseif(dir == DIR_Y) then
       print *, "Y DIRECTION"
       call write_cell_info(i,j-3,k)
       call write_cell_info(i,j-2,k)
       call write_cell_info(i,j-1,k)
       call write_cell_info(i,j+0,k)
       call write_cell_info(i,j+1,k)
       call write_cell_info(i,j+2,k)
    else
       print *, "Z DIRECTION"
       call write_cell_info(i,j,k-3)
       call write_cell_info(i,j,k-2)
       call write_cell_info(i,j,k-1)
       call write_cell_info(i,j,k+0)
       call write_cell_info(i,j,k+1)
       call write_cell_info(i,j,k+2)
    end if

    deallocate(xcent)
    deallocate(ycent)
    deallocate(zcent)


    call Driver_abortFlash( &
         "[hy_uhd_getFaceFlux]: Imaginary sound speed has obtained in " // &
         trim(solver) // " solver. " // &
         "Please try other (more diffusive) slope limiter, flux, order, cfl, etc. "//&
         "in order to increase numerical stability. LOOK AT THE LOG FILE")

  end subroutine do_error


  subroutine write_cell_info(i, j, k)
    implicit none

    integer, intent(in) :: i, j, k

    print *, ""
    print *, "CELL: ", i, j, k
    print *, "POSITION: ", xcent(i), ycent(j), zcent(k)
    print *, "DENS: ", U(DENS_VAR,i, j, k)
    print *, "PRES: ", U(PRES_VAR,i, j, k)
    print *, "GAMC: ", U(GAMC_VAR,i, j, k)

  end subroutine write_cell_info

  subroutine faceStatesEos(V,t_start,spc)

    use Eos_interface,  ONLY : Eos, Eos_getAbarZbarArraySection
    implicit none

    real, dimension(HY_VARINUMMAX), intent(INOUT) :: V
    real, intent(IN) :: t_start
    real, OPTIONAL, dimension(:), intent(IN) :: spc !dimension(HY_NSPEC)

    real, dimension(EOS_NUM) :: eosData
    integer :: interp_eosMode = MODE_DENS_PRES

    eosData(EOS_DENS) = V(HY_DENS) 
    eosData(EOS_PRES) = V(HY_PRES)
    eosData(EOS_EINT) = V(HY_EINT)
    eosData(EOS_GAMC) = V(HY_GAMC)
    ! initial guess using the previous time step temp
    eosData(EOS_TEMP) = t_start

    if (present(spc)) then
       call Eos_getAbarZbarArraySection(SPECIES_BEGIN,spc,&
            abar=eosData(EOS_ABAR),&
            zbar=eosData(EOS_ZBAR))
       call Eos(interp_eosMode,1,eosData,spc(1:NSPECIES))
    else 
       call Eos(interp_eosMode,1,eosData)
    end if

    V(HY_GAMC) = eosData(EOS_GAMC)
    V(HY_EINT) = eosData(EOS_EINT)
    V(HY_GAME) = 1.+eosData(EOS_PRES)/(eosData(EOS_DENS)*eosData(EOS_EINT)) ! game
    V(HY_PRES) = eosData(EOS_PRES)

  end subroutine faceStatesEOS

End Subroutine hy_uhd_getFaceFlux
