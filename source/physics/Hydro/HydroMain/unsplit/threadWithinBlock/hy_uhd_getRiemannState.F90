!!****if* source/physics/Hydro/HydroMain/unsplit/threadWithinBlock/hy_uhd_getRiemannState
!!
!! NAME
!!
!!  hy_uhd_getRiemannState
!!
!! SYNOPSIS
!!
!!  hy_uhd_getRiemannState( integer(IN) :: blockID,
!!                          integer(IN) :: blkLimits,
!!                          integer(IN) :: blkLimitsGC,
!!                          real(IN)    :: dt,
!!                          integer(IN) :: del(MDIM),
!!                          real(IN)    :: ogravX(:,:,:),
!!                          real(IN)    :: ogravY(:,:,:),
!!                          real(IN)    :: ogravZ(:,:,:),
!!                          real, pointer, dimension(:,:,:,:) :: scrchFaceXPtr,
!!                          real, pointer, dimension(:,:,:,:) :: scrchFaceYPtr,
!!                          real, pointer, dimension(:,:,:,:) :: scrchFaceZPtr,
!!                          real, pointer, optional, dimension(:,:,:,:,:) :: hy_SpcR,hy_SpcL,hy_SpcSig,
!!                          logical(IN), optional :: normalFieldUpdate)
!!
!! DESCRIPTION
!!
!!  This routine computes the Riemann state values at cell interfaces using 
!!  the cell centered variables and store them in the scratch arrays.
!!
!!  A 2D Cartesian configuration of a single cell is shown:
!!
!!
!!           ---------------------
!!          |          yp         |
!!          |                     |
!!          |                     |
!!          |                     |
!!          |                     |
!!          |xm      (i,j)      xp|
!!          |                     |
!!          |                     |
!!          |                     |
!!  y       |                     |
!!  |       |          ym         |
!!  |        ---------------------
!!  |______x
!!
!!
!! ARGUMENTS
!!
!!  blockID     - local block ID
!!  blkLimits   - an array that holds the lower and upper indices of the section
!!                of block without the guard cells
!!  blkLimitsGC - an array that holds the lower and upper indices of the section
!!                of block with the guard cells
!!  dt          - a current time step
!!  del         - deltas in each direction 
!!  ogravX       - gravity component in x-direction at n step
!!  ogravY       - gravity component in y-direction at n step
!!  ogravZ       - gravity component in z-direction at n step
!!  scrchFaceXPtr,scrchFaceYPtr,scrchFaceZPtr - Pointers to the scrch array (for left/right states)
!!  hy_SpcR,hy_SpcL,hy_SpcSig - Pointers for Species and mass scalar recon.
!!  normalFieldUpdate - a logical switch to choose normal magnetic fields updates only
!!                      (needed for MHD only)
!!   
!!***

!!REORDER(4):U, V0, scrchFaceXPtr,scrchFaceYPtr,scrchFaceZPtr, B[xyz]

Subroutine hy_uhd_getRiemannState(blockID,blkLimits,blkLimitsGC,dt,del,&
                                  ogravX,ogravY,ogravZ,&
                                  scrchFaceXPtr,scrchFaceYPtr,scrchFaceZPtr,&
                                  hy_SpcR,hy_SpcL,hy_SpcSig,&
                                  normalFieldUpdate)
#include "Flash.h"
  use Hydro_data,           ONLY : hy_shockDetectOn,     &
                                   hy_order,             &
                                   hy_flattening,        &
                                   hy_useGravity,        &
                                   hy_useGravHalfUpdate, &
                                   hy_fallbackLowerCFL,  &
                                   hy_use3dFullCTU,      &
                                   hy_geometry,          &
                                   hy_numXN,             &
                                   hy_cfl,               &
                                   hy_cfl_original,      &
                                   hy_cflFallbackFactor, &
                                   hy_fullSpecMsFluxHandling, &
                                   hy_useHybridOrder,    &
                                   hy_eswitch,           &
                                   hy_upwindTVD,         &
                                   hy_threadWithinBlock
#ifdef FLASH_USM_MHD
  use Hydro_data,           ONLY : hy_killDivB,   &
                                   hy_forceHydroLimit
#endif
#ifdef FLASH_UGLM_MHD
  use Hydro_data,           ONLY : hy_C_hyp, hy_C_par
#endif

  use hy_uhd_slopeLimiters, ONLY : mc
  use hy_uhd_interface,     ONLY : hy_uhd_dataReconstOnestep, &
                                   hy_uhd_shockDetect,        &
                                   hy_uhd_eigenParameters,    &
                                   hy_uhd_eigenValue,         &
                                   hy_uhd_eigenVector,        &
                                   hy_uhd_upwindTransverseFlux
  use Grid_interface,       ONLY : Grid_getBlkPtr,     &
                                   Grid_releaseBlkPtr, &
                                   Grid_getCellCoords

  implicit none

#include "constants.h"
#include "UHD.h"

  !! Arguments type declaration ------------------------------------------------------------
  integer, intent(IN)   :: blockID
  integer, intent(IN),dimension(LOW:HIGH,MDIM):: blkLimits, blkLimitsGC
  real,    intent(IN)   :: dt
  real,    intent(IN),dimension(MDIM) :: del
#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC), intent(IN) :: &
       ogravX,ogravY,ogravZ
#else
  real, dimension(blkLimitsGC(HIGH,IAXIS),blkLimitsGC(HIGH,JAXIS),blkLimitsGC(HIGH,KAXIS)), &
       intent(IN) :: ogravX,ogravY,ogravZ
#endif
  real, pointer, dimension(:,:,:,:) :: scrchFaceXPtr, scrchFaceYPtr, scrchFaceZPtr
  real, pointer, optional, dimension(:,:,:,:,:) :: hy_SpcR,hy_SpcL,hy_SpcSig
  logical, intent(IN), optional :: normalFieldUpdate
  !! ---------------------------------------------------------------------------------------

  integer :: i0,imax,j0,jmax,k0,kmax, i, j, k
  integer,dimension(MDIM) :: dataSize
  real    :: cellCfl,minCfl
  logical :: lowerCflAtBdry
  real, pointer, dimension(:,:,:,:) :: U
  integer :: dir

  real, pointer, dimension(:,:,:,:) :: Bx,By,Bz
! MHD only-------------------------------------------------------------------------------
#ifdef FLASH_USM_MHD
  logical :: normalFieldUpdateOnly
#endif
! end of ifdef FLASH_USM_MHD
! MHD only-------------------------------------------------------------------------------

#ifdef FLASH_UHD_HYDRO
  logical, save :: normalFieldUpdateOnly = .FALSE.
#endif

#ifdef FIXEDBLOCKSIZE
  real, dimension(NDIM,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: FlatCoeff,FlatTilde
  real, dimension(     GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: DivU
#else
  real, dimension(NDIM,blkLimitsGC(HIGH,IAXIS),&
                       blkLimitsGC(HIGH,JAXIS),&
                       blkLimitsGC(HIGH,KAXIS)) :: FlatCoeff, FlatTilde
  real, dimension(     blkLimitsGC(HIGH,IAXIS),&
                       blkLimitsGC(HIGH,JAXIS),&
                       blkLimitsGC(HIGH,KAXIS)) :: DivU
#endif
  real :: Sp, dv1, dp1, dp2, presL,presR,hdt

#ifdef FIXEDBLOCKSIZE  
  real, dimension(GRID_IHI_GC) :: xCenter  
  real, dimension(GRID_JHI_GC) :: yCenter  
#else  
  real, dimension(blkLimitsGC(HIGH,IAXIS)) :: xCenter  
  real, dimension(blkLimitsGC(HIGH,JAXIS)) :: yCenter  
#endif

  integer,parameter :: k2=K2D
  integer,parameter :: k3=K3D
  !! index for gravity
#ifdef GRAVITY
  integer,parameter :: kGrav = 1
#else
  integer,parameter :: kGrav = 0
#endif
  integer :: kHydro,kUSM,order
  integer :: k4,im2,ip2,jm2,jp2,km2,kp2

  ! cylindrical geometry
  integer :: velPhi, velTht, magPhi, magZ
  integer :: HY_velPhi, HY_velTht, H_magPhi, H_magZ

  real :: Rinv, geoFac, eta, enth   
  real :: sGeo_dens, sGeo_velx, sGeo_velp, sGeo_velt, sGeo_pres
  real :: sGeo_trans, sGeo_eint, sGeo_magz,sGeo_magp

  ! for 3d only --------------
  real, dimension(HY_SPEC_END) :: TransFluxXY,TransFluxYZ,TransFluxZX,&
                                  TransFluxYX,TransFluxZY,TransFluxXZ
  real :: dt2dxdy6,dt2dydz6,dt2dzdx6,hdtdx,hdtdy,hdtdz
  logical :: cons=.false.
  real    :: cs,ca,cf,as,af,uN
  !! Set transverse flux interpolation order 
  integer,parameter :: transOrder3D = 1
  real, dimension(MDIM) :: beta
  ! for 3d only --------------
  logical :: TransX_updateOnly,TransY_updateOnly,TransZ_updateOnly
  logical :: allTransUpdateOnly

  integer,dimension(4) :: lbx,ubx,lby,uby,lbz,ubz
  integer :: iDim
  real, dimension(HY_VARINUMMAX,NDIM) :: Wp, Wn
  real, dimension(HY_END_VARS) :: Vc
  real, pointer, dimension(:)   :: SigmPtr,SigcPtr,SigpPtr
  !real, allocatable,dimension(:,:,:) :: DivU

  real, allocatable,dimension(:,:,:,:,:),target :: sig
  real, allocatable,dimension(:,:,:,:,:) :: lambda
  real, allocatable,dimension(:,:,:,:,:,:) :: leftEig
  real, allocatable,dimension(:,:,:,:,:,:) :: rghtEig


#ifdef FLASH_UHD_HYDRO
  !PT: Add these variables so we can maintain a single omp parallel
  !directive for unsplit hydro and MHD simulations.  They are not used.
  logical, save :: hy_killDivB = .false.
  logical, save :: hy_forceHydroLimit = .false.
#endif
  
  
  ! GLM fluxes & state updates
#ifdef FLASH_UGLM_MHD
#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: &
       GLMxStar,GLMyStar,GLMzStar,BxStar,ByStar,BzStar
#else
  real, dimension(blkLimitsGC(HIGH,IAXIS),&
                  blkLimitsGC(HIGH,JAXIS),&
                  blkLimitsGC(HIGH,KAXIS)) :: &
       GLMxStar,GLMyStar,GLMzStar,BxStar,ByStar,BzStar
#endif
#endif /* FLASH_UGLM_MHD */




  !! index for pure Hydro
  kHydro = 1
#ifdef FLASH_USM_MHD
  if (.NOT. hy_forceHydroLimit) kHydro = 0

  normalFieldUpdateOnly = normalFieldUpdate ! the dummy argument is required to be present for MHD.
#endif
  kUSM = 1 - kHydro

  !! indices for various purposes
#if NDIM == 1
  i0   = blkLimits(LOW, IAXIS)
  imax = blkLimits(HIGH,IAXIS)
  j0   = 3
  jmax =-1
  k0   = 3
  kmax =-1
#elif NDIM == 2
  i0   = blkLimits(LOW, IAXIS)
  imax = blkLimits(HIGH,IAXIS)
  j0   = blkLimits(LOW, JAXIS)
  jmax = blkLimits(HIGH,JAXIS)
  k0   = 3
  kmax =-1
#elif NDIM == 3
  i0   = blkLimits(LOW, IAXIS)
  imax = blkLimits(HIGH,IAXIS)
  j0   = blkLimits(LOW, JAXIS)
  jmax = blkLimits(HIGH,JAXIS)
  k0   = blkLimits(LOW, KAXIS)
  kmax = blkLimits(HIGH,KAXIS)
#endif

  ! half delta t
  hdt = 0.5*dt

!!$#ifdef FLASH_UGLM_MHD
!!$  hy_C_hyp = 0.8/dt*min(del(DIR_X),del(DIR_Y))
!!$#endif

!!$  !Initialize geometric source terms
!!$  sGeo_dens=0.
!!$  sGeo_velx=0. 
!!$  sGeo_velp=0.
!!$  sGeo_pres=0.
!!$  sGeo_trans=0.
!!$  sGeo_eint=0.
!!$  sGeo_magz=0.
!!$  sGeo_magp=0.

  !! Get block pointer to UNK data
  call Grid_getBlkPtr(blockID,U,CENTER)

  !$omp parallel if (hy_threadWithinBlock) &
  !$omp default(none) &
  !$omp shared(i0,j0,k0,imax,jmax,kHydro,kUSM,kmax,&
  !$omp scrchFaceXPtr,scrchFaceYPtr,scrchFaceZPtr,&
  !$omp hy_fullSpecMsFluxHandling,normalFieldUpdateOnly,&
  !$omp blockID,blkLimitsGC,dt,del,ogravX,ogravY,ogravZ,FlatCoeff,hy_SpcR,hy_SpcL,hy_SpcSig,&
  !$omp hy_useGravity,U,Bx,By,Bz,hy_useGravHalfUpdate,leftEig,rghtEig,lambda,&
  !$omp hy_order,hy_flattening,FlatTilde,blkLimits,hdt,hy_shockDetectOn,&
  !$omp hy_use3dFullCTU,hy_geometry,hy_useHybridOrder,hy_cfl_original,&
  !$omp hy_cfl,hy_cflFallbackFactor,hy_fallbackLowerCFL,minCfl,&
  !$omp xCenter,yCenter,DivU,sig,hy_forceHydroLimit,hy_killdivb,datasize) &
  !$omp private(i,j,k,lbx,ubx,lby,uby,lbz,ubz,iDim, un,cf,Vc,order,Wn,Wp,sigmptr,sigcptr,sigpptr,&
  !$omp k4,im2,ip2,jm2,jp2,km2,kp2,&
  !$omp dp1,dp2,dv1,presL,presR,Sp,magPhi,magZ,sGeo_magp,sGeo_magz,&
  !$omp TransX_updateOnly,TransY_updateOnly,TransZ_updateOnly,allTransUpdateOnly,&
  !$omp cellCfl,lowerCflAtBdry,&
  !$omp dt2dxdy6,dt2dydz6,dt2dzdx6,TransFluxXY,TransFluxYZ,TransFluxZX,&
  !$omp TransFluxYX,TransFluxZY,TransFluxXZ,Rinv,velPhi,velTht,&
  !$omp sGeo_dens,sGeo_velx,sGeo_velp,sGeo_velt,sGeo_pres,sGeo_trans,sGeo_eint,&
  !$omp geoFac,cs,eta,enth,dir,HY_velPhi,HY_velTht,h_magphi,h_magz)


  ! MHD only-------------------------------------------------------------------------------
#if defined(FLASH_USM_MHD) && (NFACE_VARS > 0) && (NDIM > 1)
  if (hy_order > 1) then
     !$omp sections      ! assume we do not use FL_NON_PERMANENT_GUARDCELLS
     !$omp section
     call Grid_getBlkPtr(blockID,Bx,FACEX)
     !$omp section
     call Grid_getBlkPtr(blockID,By,FACEY)
     !$omp section
     if (NDIM == 3) call Grid_getBlkPtr(blockID,Bz,FACEZ)
     !$omp end sections
  endif
#endif /* endif of if defined(FLASH_USM_MHD) && NFACE_VARS > 0 && NDIM > 1 */
  ! MHD only-------------------------------------------------------------------------------



  if (.NOT. normalFieldUpdateOnly) then

     !$omp sections
     !$omp section
     if (hy_geometry /= CARTESIAN) then
        ! Grab cell x-coords for this block  
        call Grid_getCellCoords(IAXIS,blockID, CENTER,.true.,xCenter, blkLimitsGC(HIGH,IAXIS))
     endif
#if NDIM > 1
     !$omp section
     if (hy_geometry == SPHERICAL) &
             call Grid_getCellCoords(JAXIS,blockID, CENTER,.true.,yCenter, blkLimitsGC(HIGH,JAXIS))
#endif
     !$omp section
     dataSize(1:MDIM)=blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM)+1
     allocate( sig(HY_VARINUMMAX,NDIM,           dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
     allocate( lambda(HY_WAVENUM,NDIM,           dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
     allocate(leftEig(HY_VARINUM,HY_WAVENUM,NDIM,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
     allocate(rghtEig(HY_VARINUM,HY_WAVENUM,NDIM,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
     !$omp end sections
  end if

  !! -----------------------------------------------------------------------!
  !! (1) Compute eigen structures once and for all -------------------------!
  !! (2) Compute divergence of velocity field      -------------------------!
  !! -----------------------------------------------------------------------!
  !if (hy_upwindTVD) kHydro = -1

  if (hy_useHybridOrder .AND. .NOT. normalFieldUpdateOnly) then
     k4 = hy_order - 1             !cf. hy_uhd_dataReconstOnestep
     if (k4 > 2) k4 = 2 !(i.e., assume order = 3)
     k4 = min(NGUARD-2-kUSM,k4)
     im2=max(blkLimitsGC(LOW,IAXIS)+1 ,i0  -2+kHydro-k4)
     ip2=min(blkLimitsGC(HIGH,IAXIS)-1,imax+2-kHydro+k4)
#if NDIM > 1
     jm2=max(blkLimitsGC(LOW,JAXIS)+1 ,j0  -2+kHydro-k4)
     jp2=min(blkLimitsGC(HIGH,JAXIS)-1,jmax+2-kHydro+k4)
#else
     jm2=1; jp2=1
#endif
#if NDIM > 2
     km2=max(blkLimitsGC(LOW,KAXIS)+1 ,k0  -2+kHydro-k4)
     kp2=min(blkLimitsGC(HIGH,KAXIS)-1,kmax+2-kHydro+k4)
#else
     km2=1; kp2=1
#endif
     !$omp do schedule(static)
     do k=km2,kp2
        do j=jm2,jp2
           do i=im2,ip2
              !! We used to compute undivided divergence of velocity fields and (magneto)sonic speed
              !! and store local (magneto)sonic speeds for hybrid order. Now the latter are
              !! computed elsewhere.

              DivU(i,j,k) = U(VELX_VAR,i+1,j,k)-U(VELX_VAR,i-1,j,k)
              if (NDIM > 1) then
                 DivU(i,j,k) = DivU(i,j,k) &
                      +U(VELY_VAR,i,j+1,k)-U(VELY_VAR,i,j-1,k)
                 if (NDIM > 2) then
                    DivU(i,j,k) = DivU(i,j,k) &
                         +U(VELZ_VAR,i,j,k+1)-U(VELZ_VAR,i,j,k-1)
                 endif
              endif
              DivU(i,j,k) = 0.5*DivU(i,j,k)

           enddo ! do i-loop
        enddo ! do j-loop
     enddo ! do k-loop
     !$omp end do
  endif !end of hybridOrder


  !! -----------------------------------------------------------------------!
  !! (3) Flattening begins here --------------------------------------------!
  !! -----------------------------------------------------------------------!
  if (hy_flattening .AND. .NOT. normalFieldUpdateOnly) then

     ! Initialize with zero
     !$omp workshare     
     FlatTilde(:,:,:,:) = 0.
     FlatCoeff(:,:,:,:) = 0.
     !$omp end workshare

     ! Flat tilde
     !$omp do schedule(static)
     do k=k0-2-min(NGUARD-4,kUSM)*k3,kmax+2+min(NGUARD-4,kUSM)*k3
        do j=j0-2-min(NGUARD-4,kUSM)*k2,jmax+2+min(NGUARD-4,kUSM)*k2
           do i=i0-2-min(NGUARD-4,kUSM),imax+2+min(NGUARD-4,kUSM)
              do dir=1,NDIM

                 select case (dir)
                 case (DIR_X)
                    dp1   = (U(PRES_VAR,i+1,j,k)-U(PRES_VAR,i-1,j,k))
                    dp2   = (U(PRES_VAR,i+2,j,k)-U(PRES_VAR,i-2,j,k))
                    dv1   =  U(VELX_VAR,i+1,j,k)-U(VELX_VAR,i-1,j,k)
                    presL =  U(PRES_VAR,i-1,j,k)
                    presR =  U(PRES_VAR,i+1,j,k)
#if NDIM > 1
                 case (DIR_Y)
                    dp1   = (U(PRES_VAR,i,j+1,k)-U(PRES_VAR,i,j-1,k))
                    dp2   = (U(PRES_VAR,i,j+2,k)-U(PRES_VAR,i,j-2,k))
                    dv1   =  U(VELY_VAR,i,j+1,k)-U(VELY_VAR,i,j-1,k)
                    presL =  U(PRES_VAR,i,j-1,k)
                    presR =  U(PRES_VAR,i,j+1,k)
#if NDIM > 2
                 case (DIR_Z)
                    dp1   = (U(PRES_VAR,i,j,k+1)-U(PRES_VAR,i,j,k-1))
                    dp2   = (U(PRES_VAR,i,j,k+2)-U(PRES_VAR,i,j,k-2))
                    dv1   =  U(VELZ_VAR,i,j,k+1)-U(VELZ_VAR,i,j,k-1)
                    presL =  U(PRES_VAR,i,j,k-1)
                    presR =  U(PRES_VAR,i,j,k+1)
#endif
#endif
                 end select

                 if (abs(dp2) > 1.e-15) then
                    Sp = dp1/dp2 - 0.75
                 else
                    Sp = 0.
                 endif

                 FlatTilde(dir,i,j,k) = max(0.0, min(1.0,10.0*Sp))
                 if ((abs(dp1)/min(presL,presR) < 1./3.) .or. dv1 > 0. ) then
                    FlatTilde(dir,i,j,k) = 0.
                 endif

              enddo
           enddo
        enddo
     enddo
     !$omp end do

     ! Flat coefficient
     !$omp do schedule(static)
     do k=k0-2+kHydro*k3,kmax+2-kHydro*k3
        do j=j0-2+kHydro*k2,jmax+2-kHydro*k2
           do i=i0-2+kHydro,imax+2-kHydro
              do dir=1,NDIM

                 select case (dir)
                 case (DIR_X)
                    dp1   = (U(PRES_VAR,i+1,j,k)-U(PRES_VAR,i-1,j,k))

                    if ( dp1 < 0.0 ) then
                       FlatCoeff(dir,i,j,k) = max(FlatTilde(dir,i,j,k),FlatTilde(dir,i+1,j,k))
                    elseif (dp1 == 0.) then
                       FlatCoeff(dir,i,j,k) = FlatTilde(dir,i,j,k)
                    else
                       FlatCoeff(dir,i,j,k) = max(FlatTilde(dir,i,j,k),FlatTilde(dir,i-1,j,k))
                    endif
#if NDIM > 1
                 case (DIR_Y)
                    dp1   = (U(PRES_VAR,i,j+1,k)-U(PRES_VAR,i,j-1,k))

                    if ( dp1 < 0.0 ) then
                       FlatCoeff(dir,i,j,k) = max(FlatTilde(dir,i,j,k),FlatTilde(dir,i,j+1,k))
                    elseif (dp1 == 0.) then
                       FlatCoeff(dir,i,j,k) = FlatTilde(dir,i,j,k)
                    else
                       FlatCoeff(dir,i,j,k) = max(FlatTilde(dir,i,j,k),FlatTilde(dir,i,j-1,k))
                    endif
#if NDIM > 2
                 case (DIR_Z)
                    dp1   = (U(PRES_VAR,i,j,k+1)-U(PRES_VAR,i,j,k-1))

                    if ( dp1 < 0.0 ) then
                       FlatCoeff(dir,i,j,k) = max(FlatTilde(dir,i,j,k),FlatTilde(dir,i,j,k+1))
                    elseif (dp1 == 0.) then
                       FlatCoeff(dir,i,j,k) = FlatTilde(dir,i,j,k)
                    else
                       FlatCoeff(dir,i,j,k) = max(FlatTilde(dir,i,j,k),FlatTilde(dir,i,j,k-1))
                    endif
#endif
#endif
                 end select
              enddo
           enddo
        enddo
     enddo
     !$omp end do
  endif

  !! -----------------------------------------------------------------------!
  !! (4) Start calculating Riemann states ----------------------------------!
  !! -----------------------------------------------------------------------!
  !! Compute Riemann states at each cell
  if (.not. normalFieldUpdateOnly) then

#ifdef CFL_VAR
     minCfl = hy_cfl_original
#else
     minCfl = hy_cfl
#endif

     !$omp do reduction(min:minCfl) schedule(static)
     do k=k0-2-(k3*kUSM-kHydro)*k3,kmax+2+(k3*kUSM-kHydro)*k3

        do j=j0-2-(k3*kUSM-kHydro)*k2,jmax+2+(k3*kUSM-kHydro)*k2
           do i=i0-2-(k3*kUSM-kHydro),imax+2+(k3*kUSM-kHydro)
           ! Extra stencil is needed for 3D to correctly calculate transverse fluxes 
           !(i.e., cross derivatives in x,y, & z)

              !! save the cell center values for later use
              Vc(HY_DENS:HY_END_VARS-kGrav) = &
                      (/U(DENS_VAR,i,j,k)&
                       ,U(VELX_VAR:VELZ_VAR,i,j,k)&
                       ,U(PRES_VAR,i,j,k)&
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
                       ,U(MAGX_VAR:MAGZ_VAR,i,j,k)&
#endif
#ifdef FLASH_UGLM_MHD
                       ,U(GLMP_VAR,i,j,k) &
#endif
                       ,U(GAMC_VAR,i,j,k) &
                       ,U(GAME_VAR,i,j,k) &
                       ,U(EINT_VAR,i,j,k) &
#ifdef FLASH_UHD_3T
                       ,U(EELE_VAR,i,j,k) &
                       ,U(EION_VAR,i,j,k) &
                       ,U(ERAD_VAR,i,j,k) &
#endif
                       /)

              order = hy_order
              lowerCflAtBdry = .FALSE.
#ifdef BDRY_VAR
              !! Reduce order in fluid cell near solid boundary if defined:
              !! Reduce order of spatial reconstruction depending on the distance to the solid boundary
              if (order > 2) then
                 im2=max(blkLimitsGC(LOW,IAXIS),i-2); ip2=min(blkLimitsGC(HIGH,IAXIS),i+2)
#if NDIM > 1
                 jm2=max(blkLimitsGC(LOW,JAXIS),j-2); jp2=min(blkLimitsGC(HIGH,JAXIS),j+2)
#else
                 jm2 = 1; jp2=1
#endif
#if NDIM > 2
                 km2=max(blkLimitsGC(LOW,KAXIS),k-2); kp2=min(blkLimitsGC(HIGH,KAXIS),k+2)
#else
                 km2 = 1; kp2=1
#endif
                 if (maxval(U(BDRY_VAR,im2:ip2,jm2:jp2,km2:kp2)) .LE. 0.) then !everyone is fluid
                    order = 3
                 else
                    order = 2
                 endif
              endif
              if ((U(BDRY_VAR,i,j,k)*U(BDRY_VAR,i-1, j,   k   ) < 0.0) .or. &
                  (U(BDRY_VAR,i,j,k)*U(BDRY_VAR,i+1, j,   k   ) < 0.0) .or. &
                  (U(BDRY_VAR,i,j,k)*U(BDRY_VAR,i,   j-k2,k   ) < 0.0) .or. &
                  (U(BDRY_VAR,i,j,k)*U(BDRY_VAR,i,   j+k2,k   ) < 0.0) .or. &
                  (U(BDRY_VAR,i,j,k)*U(BDRY_VAR,i,   j,   k-k3) < 0.0) .or. &
                  (U(BDRY_VAR,i,j,k)*U(BDRY_VAR,i,   j,   k+k3) < 0.0)) then
                 order = 1
                 lowerCflAtBdry = .TRUE.
              endif
#if NDIM > 2
              ! Addtional 3D test whether a solid cell is so close
              ! that we cannot do all the proper transverse
              ! computations for this cell.
              if (hy_use3dFullCTU) then
                 if (maxval(U(BDRY_VAR,i-1:i+1,j-1:j+1,k-1:k+1)) > 0.) then
                    if (maxval(U(BDRY_VAR,i      ,j-1:j+1,k-1:k+1)) > 0.) lowerCflAtBdry = .TRUE.
                    if (maxval(U(BDRY_VAR,i-1:i+1,j      ,k-1:k+1)) > 0.) lowerCflAtBdry = .TRUE.
                    if (maxval(U(BDRY_VAR,i-1:i+1,j-1:j+1,k      )) > 0.) lowerCflAtBdry = .TRUE.
                 end if
              end if
#endif
#endif

#ifdef CFL_VAR
              cellCfl = U(CFL_VAR,i,j,k)
#else
              cellCfl = hy_cfl
#endif
              if (lowerCflAtBdry) cellCfl = min(cellCfl, hy_cflFallbackFactor / real(NDIM))

              !! Flag for tranverse update 
              TransX_updateOnly = .false.
              TransY_updateOnly = .false.
              TransZ_updateOnly = .false.

              if (i > blkLimits(HIGH,IAXIS)+1 .or. &
                   i < blkLimits(LOW, IAXIS)-1) then
                 TransX_updateOnly = .true.
              endif
              if (i > blkLimits(HIGH,IAXIS)+2*kUSM .or. &
                   i < blkLimits(LOW, IAXIS)-2*kUSM) then
                 TransY_updateOnly = .true.
                 TransZ_updateOnly = .true.
              endif
#if NDIM > 1
              if (j > blkLimits(HIGH,JAXIS)+1 .or. &
                   j < blkLimits(LOW, JAXIS)-1) then
                 TransY_updateOnly = .true.
              endif
              if (j > blkLimits(HIGH,JAXIS)+2*kUSM .or. &
                   j < blkLimits(LOW, JAXIS)-2*kUSM) then
                 TransX_updateOnly = .true.
                 TransZ_updateOnly = .true.
              endif

#if NDIM > 2
              if (k > blkLimits(HIGH,KAXIS)+1 .or. & 
                   k < blkLimits(LOW, KAXIS)-1) then
                 TransZ_updateOnly = .true.
              endif
              if (k > blkLimits(HIGH,KAXIS)+2*kUSM .or. & 
                   k < blkLimits(LOW, KAXIS)-2*kUSM) then
                 TransX_updateOnly = .true.
                 TransY_updateOnly = .true.
              endif
#endif
#endif
              allTransUpdateOnly = TransX_updateOnly
              if (NDIM > 1) allTransUpdateOnly = allTransUpdateOnly .AND. TransY_updateOnly
              if (NDIM > 2) allTransUpdateOnly = allTransUpdateOnly .AND. TransZ_updateOnly

              if (order == 1) then
                 !! DEV: THE FIRST ORDER SHOULD GO INTO THE DATA RECONSTRUCT ONE STEP TO GET
                 !! TRANSVERSE FLUXES
                 if (.NOT. allTransUpdateOnly) then
                    do iDim = 1,NDIM
                       call fallbackToFirstOrder(iDim,Wn(:,iDim),Wp(:,iDim),Vc,hy_SpcL,hy_SpcR,U,i,j,k)
                    enddo
                 end if

              else ! for high-order schemes


                 !! Left and right Riemann state reconstructions
                 if (hy_fullSpecMsFluxHandling .AND. hy_numXN > 0 &
                      .AND. present(hy_spcR)) then
                    call hy_uhd_dataReconstOnestep&
                      (blockID,blkLimitsGC,    &
                       order,i,j,k,dt,del,     &
                       ogravX,ogravY,ogravZ,   &
                       DivU,FlatCoeff,         &
                       TransX_updateOnly,      &
                       TransY_updateOnly,      &
                       TransZ_updateOnly,      &
                       Wp,Wn,                  &
                       sig    (  1,1,i,j,k),   &
                       lambda (  1,1,i,j,k),   &
                       leftEig(1,1,1,i,j,k),   &
                       rghtEig(1,1,1,i,j,k),   &
                       cellCfl, &
                       hy_SpcR,hy_SpcL,hy_SpcSig)
                 else
                    call hy_uhd_dataReconstOnestep&
                      (blockID,blkLimitsGC,    &
                       order,i,j,k,dt,del,     &
                       ogravX,ogravY,ogravZ,   &
                       DivU,FlatCoeff,         &
                       TransX_updateOnly,      &
                       TransY_updateOnly,      &
                       TransZ_updateOnly,      &
                       Wp,Wn,                  &
                       sig    (  1,1,i,j,k),   &
                       lambda (  1,1,i,j,k),   &
                       leftEig(1,1,1,i,j,k),   &
                       rghtEig(1,1,1,i,j,k),   & 
                       cellCfl)
                 endif ! if(hy_fullSpecMsFluxHandling ...

              endif ! end of high-order reconstruction schemes


              if (hy_geometry /= CARTESIAN) then
                 !! **************************************************************
                 !! Add geometric source terms in left and Right States          *
                 !! **************************************************************


                 !Initialize geometric source terms
                 sGeo_dens = 0.
                 sGeo_velx = 0. 
                 sGeo_velp = 0.
                 sGeo_pres = 0.
                 sGeo_eint = 0.
                 sGeo_magz = 0.
                 sGeo_magp = 0.
                 sGeo_trans= 0.

                 Rinv = 1./xCenter(i)
                 select case (hy_geometry)
                 case (CYLINDRICAL)
                    velPhi    = VELZ_VAR
                    HY_velPhi = HY_VELZ
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
                    magPhi    = MAGZ_VAR
                    magZ      = MAGY_VAR
                    H_magPhi  = HY_MAGZ
                    H_magZ    = HY_MAGY
#endif
                    geoFac    = Rinv
                 case (POLAR)
                    velPhi    = VELY_VAR
                    HY_velPhi = HY_VELY
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
                    magPhi    = MAGY_VAR
                    magZ      = MAGZ_VAR
                    H_magPhi  = HY_MAGY
                    H_magZ    = HY_MAGZ
#endif
                    geoFac    = Rinv
                 case (SPHERICAL)
                    velPhi    = VELZ_VAR
                    velTht    = VELY_VAR
                    HY_velPhi = HY_VELZ
                    HY_velTht = HY_VELY
                    geoFac    = 2.*Rinv
                 end select

                 cs  = sqrt(U(GAMC_VAR,i,j,k)*U(PRES_VAR,i,j,k)/U(DENS_VAR,i,j,k))
                 eta = (abs(U(VELX_VAR,i,j,k)) + cs) * dt/del(DIR_X)
                 eta = (1.-eta) / (cs*dt*abs(geoFac))
                 eta = min(1.,eta)
                 !! comment this line not to use the axis hack
!!$                 geoFac = eta * geoFac
                 !! end of the axis hack
                 enth = U(EINT_VAR,i,j,k) + U(PRES_VAR,i,j,k)/U(DENS_VAR,i,j,k)

                 !! right/left state source terms
                 sGeo_dens = -U(DENS_VAR,i,j,k) * U(VELX_VAR,i,j,k) * geoFac !src[DN]
#if NDIM > 1
                 if (hy_geometry == SPHERICAL) then
                      sGeo_dens = sGeo_dens &
                      -U(DENS_VAR,i,j,k)*U(velTht,i,j,k)*cos(yCenter(j))/sin(yCenter(j)) * 0.5*geoFac
                      sGeo_velt = (U(velPhi,i,j,k)**2) * cos(yCenter(j))/sin(yCenter(j)) * 0.5*geoFac
                   end if
#endif
                 sGeo_velx = (U(velPhi,i,j,k)**2) * geoFac                   !src[VR]  
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
                 sGeo_velx = sGeo_velx - (U(magPhi,i,j,k)**2) * geoFac / U(DENS_VAR,i,j,k)
                 if (hy_geometry == SPHERICAL) &
                      sGeo_velx = sGeo_velx + U(velTht,i,j,k)**2 * geoFac
#endif                  
                 sGeo_velp = -U(velPhi,i,j,k) * U(VELX_VAR,i,j,k) * geoFac !src[Vphi]
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
                 sGeo_velp = sGeo_velp + U(magPhi,i,j,k) * U(MAGX_VAR,i,j,k) * geoFac / U(DENS_VAR,i,j,k)
#endif                  
                 sGeo_pres = sGeo_dens * cs**2                              !src[PR]
                 sGeo_eint = (sGeo_dens*enth)/U(DENS_VAR,i,j,k)  
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
                 sGeo_magp = - U(velPhi,  i,j,k) * U(MAGX_VAR,i,j,k) * geoFac 
                 sGeo_magz = - U(VELX_VAR,i,j,k) * U(magZ,    i,j,k) * geoFac 
#endif


                 !! Add sources terms for n+1/2 Left state
                 Wn(HY_DENS,  DIR_X) = Wn(HY_DENS,  DIR_X) + hdt * sGeo_dens
                 Wn(HY_VELX,  DIR_X) = Wn(HY_VELX,  DIR_X) + hdt * sGeo_velx
                 Wn(HY_velPhi,DIR_X) = Wn(HY_velPhi,DIR_X) + hdt * sGeo_velp
                 Wn(HY_PRES,  DIR_X) = Wn(HY_PRES,  DIR_X) + hdt * sGeo_pres
                 Wn(HY_EINT,  DIR_X) = Wn(HY_EINT,  DIR_X) + hdt * sGeo_eint
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
                 Wn(H_magPhi, DIR_X) = Wn(H_magPhi, DIR_X) + hdt * sGeo_magp
                 Wn(H_magZ,   DIR_X) = Wn(H_magZ,   DIR_X) + hdt * sGeo_magz
#endif                 
                 !! Add source terms for n+1/2 Right state
                 Wp(HY_DENS,  DIR_X) = Wp(HY_DENS,  DIR_X) + hdt * sGeo_dens
                 Wp(HY_VELX,  DIR_X) = Wp(HY_VELX,  DIR_X) + hdt * sGeo_velx
                 Wp(HY_velPhi,DIR_X) = Wp(HY_velPhi,DIR_X) + hdt * sGeo_velp
                 Wp(HY_PRES,  DIR_X) = Wp(HY_PRES,  DIR_X) + hdt * sGeo_pres
                 Wp(HY_EINT,  DIR_X) = Wp(HY_EINT,  DIR_X) + hdt * sGeo_eint
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
                 Wp(H_magPhi, DIR_X) = Wp(H_magPhi, DIR_X) + hdt * sGeo_magp
                 Wp(H_magZ,   DIR_X) = Wp(H_magZ,   DIR_X) + hdt * sGeo_magz
#endif                 


                 if (xCenter(i) - 0.5*del(DIR_X) == 0.) then 
                    ! the velocity should be zero at r=0.
                    Wn(HY_VELX,  DIR_X) = 0.0
                    Wn(HY_velPhi,DIR_X) = 0.0
                    if (hy_geometry == SPHERICAL) Wn(HY_velTht,DIR_X) = 0.0
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
                    Wn(HY_MAGX,  DIR_X) = 0.0
                    Wn(H_magPhi, DIR_X) = 0.0
#endif                 
                 elseif (xCenter(i) + 0.5*del(DIR_X) == 0.) then
                    Wp(HY_VELX,  DIR_X) = 0.0
                    Wp(HY_velPhi,DIR_X) = 0.0
                    if (hy_geometry == SPHERICAL) Wp(HY_velTht,DIR_X) = 0.0
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
                    Wp(HY_MAGX,  DIR_X) = 0.0
                    Wp(H_magPhi, DIR_X) = 0.0
#endif                 
                 endif
                 !! Calculate R-momentum geometric source term for transverse fluxes
                 !! We will use the cell-centered states at t^n
                 sGeo_trans = (U(velPhi,i,j,k)**2) * geoFac
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
                 sGeo_trans = sGeo_trans - (U(magPhi,i,j,k)**2) * geoFac/U(DENS_VAR,i,j,k)
#endif
                 if (hy_geometry == SPHERICAL) &
                      sGeo_trans = sGeo_trans + U(velTht,i,j,k)**2 * geoFac
#if NDIM > 1
                 Wn(HY_VELX,DIR_Y) = Wn(HY_VELX,DIR_Y) + hdt * sGeo_trans
                 Wp(HY_VELX,DIR_Y) = Wp(HY_VELX,DIR_Y) + hdt * sGeo_trans
#if NDIM == 3
                 Wn(HY_VELX,DIR_Z) = Wn(HY_VELX,DIR_Z) + hdt * sGeo_trans
                 Wp(HY_VELX,DIR_Z) = Wp(HY_VELX,DIR_Z) + hdt * sGeo_trans
#endif
#endif
                 !! Check positivity of density and pressure
                 if (Wn(HY_DENS,DIR_X) < 0. .or. Wp(HY_DENS,DIR_X) < 0. .or. &
                     Wn(HY_PRES,DIR_X) < 0. .or. Wp(HY_PRES,DIR_X) < 0.) then
                    if(.NOT.TransX_updateOnly) &
                         call fallbackToFirstOrder(DIR_X,Wn(:,DIR_X),Wp(:,DIR_X),Vc,hy_SpcL,hy_SpcR,U,i,j,k)
                    cellCfl = min(cellCfl, hy_cflFallbackFactor / real(NDIM))
                 end if
#if NDIM > 1
                 if (Wn(HY_DENS,DIR_Y) < 0. .or. Wp(HY_DENS,DIR_Y) < 0. .or. &
                     Wn(HY_PRES,DIR_Y) < 0. .or. Wp(HY_PRES,DIR_Y) < 0.) then
                    if(.NOT.TransY_updateOnly) &
                         call fallbackToFirstOrder(DIR_Y,Wn(:,DIR_Y),Wp(:,DIR_Y),Vc,hy_SpcL,hy_SpcR,U,i,j,k)
                    cellCfl = min(cellCfl, hy_cflFallbackFactor / real(NDIM))
                 end if
#if NDIM ==3
                 if (Wn(HY_DENS,DIR_Z) < 0. .or. Wp(HY_DENS,DIR_Z) < 0. .or. &
                     Wn(HY_PRES,DIR_Z) < 0. .or. Wp(HY_PRES,DIR_Z) < 0.) then
                    if(.NOT.TransZ_updateOnly) &
                         call fallbackToFirstOrder(DIR_Z,Wn(:,DIR_Z),Wp(:,DIR_Z),Vc,hy_SpcL,hy_SpcR,U,i,j,k)
                    cellCfl = min(cellCfl, hy_cflFallbackFactor / real(NDIM))
                 end if
#endif
#endif
              endif !end if of if (hy_geometry .ne. CARTESIAN)

#ifdef CFL_VAR
              U(CFL_VAR,i,j,k) = cellCfl
#endif
              minCfl = min(minCfl,cellCfl)


#ifdef GRAVITY
              if (hy_useGravity .and. hy_useGravHalfUpdate)then
                 Wp(HY_VELX:HY_VELZ,DIR_X)=Wp(HY_VELX:HY_VELZ,DIR_X)&
                      +hdt*(/Wp(HY_GRAV,DIR_X),ogravY(i,j,k),ogravZ(i,j,k)/)

                 Wn(HY_VELX:HY_VELZ,DIR_X)=Wn(HY_VELX:HY_VELZ,DIR_X)&
                      +hdt*(/Wn(HY_GRAV,DIR_X),ogravY(i,j,k),ogravZ(i,j,k)/)
              endif
#endif  
              !! Store Riemann states to scratch arrays

              if(.NOT.TransX_updateOnly) scrchFaceXPtr(HY_P01_FACEXPTR_VAR:HY_P01_FACEXPTR_VAR+HY_SCRATCH_NUM-1,i,j,k)&
                   = Wp(HY_DENS:HY_END_VARS-kGrav,DIR_X)

              if(.NOT.TransX_updateOnly) scrchFaceXPtr(HY_N01_FACEXPTR_VAR:HY_N01_FACEXPTR_VAR+HY_SCRATCH_NUM-1,i,j,k)&
                   = Wn(HY_DENS:HY_END_VARS-kGrav,DIR_X)

#if NDIM >= 2
#ifdef GRAVITY
              if (hy_useGravity .and. hy_useGravHalfUpdate)then
                 Wp(HY_VELX:HY_VELZ,DIR_Y)=Wp(HY_VELX:HY_VELZ,DIR_Y)&
                      +hdt*(/ogravX(i,j,k),Wp(HY_GRAV,DIR_Y),ogravZ(i,j,k)/)

                 Wn(HY_VELX:HY_VELZ,DIR_Y)=Wn(HY_VELX:HY_VELZ,DIR_Y)&
                      +hdt*(/ogravX(i,j,k),Wn(HY_GRAV,DIR_Y),ogravZ(i,j,k)/)
              endif
#endif
              !! Store Riemann states to scratch arrays
              if(.NOT.TransY_updateOnly) scrchFaceYPtr(HY_P01_FACEYPTR_VAR:HY_P01_FACEYPTR_VAR+HY_SCRATCH_NUM-1,i,j,k)&
                   = Wp(HY_DENS:HY_END_VARS-kGrav,DIR_Y)

              if(.NOT.TransY_updateOnly) scrchFaceYPtr(HY_N01_FACEYPTR_VAR:HY_N01_FACEYPTR_VAR+HY_SCRATCH_NUM-1,i,j,k)&
                   = Wn(HY_DENS:HY_END_VARS-kGrav,DIR_Y)
#if NDIM == 3
#ifdef GRAVITY
              if (hy_useGravity .and. hy_useGravHalfUpdate)then
                 Wp(HY_VELX:HY_VELZ,DIR_Z)=Wp(HY_VELX:HY_VELZ,DIR_Z)&
                      +hdt*(/ogravX(i,j,k),ogravY(i,j,k),Wp(HY_GRAV,DIR_Z)/)

                 Wn(HY_VELX:HY_VELZ,DIR_Z)=Wn(HY_VELX:HY_VELZ,DIR_Z)&
                      +hdt*(/ogravX(i,j,k),ogravY(i,j,k),Wn(HY_GRAV,DIR_Z)/)
              endif
#endif
              !! Store Riemann states to scratch arrays
              if(.NOT.TransZ_updateOnly) scrchFaceZPtr(HY_P01_FACEZPTR_VAR:HY_P01_FACEZPTR_VAR+HY_SCRATCH_NUM-1,i,j,k)&
                   = Wp(HY_DENS:HY_END_VARS-kGrav,DIR_Z)

              if(.NOT.TransZ_updateOnly) scrchFaceZPtr(HY_N01_FACEZPTR_VAR:HY_N01_FACEZPTR_VAR+HY_SCRATCH_NUM-1,i,j,k)&
                   = Wn(HY_DENS:HY_END_VARS-kGrav,DIR_Z)
#endif
#endif

#if NFACE_VARS > 0
#if NDIM > 1
              if (hy_order > 1 .and. hy_killdivB .and. (.not. hy_forceHydroLimit)) then
                 if(.NOT.TransX_updateOnly) scrchFaceXPtr(HY_P06_FACEXPTR_VAR,i,j,k)= Bx(MAG_FACE_VAR,i+1, j,   k  )
                 if(.NOT.TransX_updateOnly) scrchFaceXPtr(HY_N06_FACEXPTR_VAR,i,j,k)= Bx(MAG_FACE_VAR,i,   j,   k  )

                 if(.NOT.TransY_updateOnly) scrchFaceYPtr(HY_P07_FACEYPTR_VAR,i,j,k)= By(MAG_FACE_VAR,i,   j+1, k  )
                 if(.NOT.TransY_updateOnly) scrchFaceYPtr(HY_N07_FACEYPTR_VAR,i,j,k)= By(MAG_FACE_VAR,i,   j,   k  )
#if NDIM == 3
                 if(.NOT.TransZ_updateOnly) scrchFaceZPtr(HY_P08_FACEZPTR_VAR,i,j,k)= Bz(MAG_FACE_VAR,i,   j,   k+1)
                 if(.NOT.TransZ_updateOnly) scrchFaceZPtr(HY_N08_FACEZPTR_VAR,i,j,k)= Bz(MAG_FACE_VAR,i,   j,   k  )
#endif
              endif
#endif
#endif
           enddo
        enddo
     enddo
  !$omp end do 

! MHD only-------------------------------------------------------------------------------
#ifdef FLASH_USM_MHD
  else ! else of if (.not. normalFieldUpdateOnly) then
       !! Updates of magnetic fields in normal direction
#if NFACE_VARS > 0
#if NDIM > 1
     if (hy_order > 1 .and. hy_killdivB .and. (.not. hy_forceHydroLimit)) then
        lbx = lbound(scrchFaceXPtr); ubx = ubound(scrchFaceXPtr)
        lby = lbound(scrchFaceYPtr); uby = ubound(scrchFaceYPtr)
        if (NDIM == 3) then
           lbz = lbound(scrchFaceZPtr); ubz = ubound(scrchFaceZPtr)
        end if
        !$omp workshare
        scrchFaceXPtr(HY_P06_FACEXPTR_VAR,:,:,:) = Bx(MAGI_FACE_VAR, lbx(2)+1:ubx(2)+1, lbx(3)  :ubx(3),   lbx(4):ubx(4))
        scrchFaceXPtr(HY_N06_FACEXPTR_VAR,:,:,:) = Bx(MAGI_FACE_VAR, lbx(2)  :ubx(2),   lbx(3)  :ubx(3),   lbx(4):ubx(4))
        scrchFaceYPtr(HY_P07_FACEYPTR_VAR,:,:,:) = By(MAGI_FACE_VAR, lby(2)  :uby(2),   lby(3)+1:uby(3)+1, lby(4):uby(4))
        scrchFaceYPtr(HY_N07_FACEYPTR_VAR,:,:,:) = By(MAGI_FACE_VAR, lby(2)  :uby(2),   lby(3)  :uby(3),   lby(4):uby(4))
#if NDIM == 3
        scrchFaceZPtr(HY_P08_FACEZPTR_VAR,:,:,:) = Bz(MAGI_FACE_VAR, lbz(2) :ubz(2),  lbz(3) :ubz(3),  lbz(4)+1:ubz(4)+1)
        scrchFaceZPtr(HY_N08_FACEZPTR_VAR,:,:,:) = Bz(MAGI_FACE_VAR, lbz(2) :ubz(2),  lbz(3) :ubz(3),  lbz(4)  :ubz(4))
#endif
        !$omp end workshare 
     endif
#endif
#endif
#endif /* end of ifdef FLASH_USM_MHD */
! MHD only-------------------------------------------------------------------------------
  endif ! end of if (.not. normalFieldUpdateOnly) then

#ifndef CFL_VAR
#if NDIM > 1
  if (.not. normalFieldUpdateOnly) then
#ifdef BDRY_VAR
     !$omp critical (hy_crit_update_hy_cfl)
     hy_cfl = minCfl
     !$omp end critical (hy_crit_update_hy_cfl)
#else
     if (hy_useHybridOrder .OR. hy_fallbackLowerCFL) then
        !$omp critical (hy_crit_update_hy_cfl)
        hy_cfl = minCfl
        !$omp end critical (hy_crit_update_hy_cfl)
     end if
#endif
  end if
#endif
#endif


  !!---------------------------------------------------------------!
  !! (5) Transverse correction terms for 3D -----------------------!
  !!---------------------------------------------------------------!
#if NDIM == 3
  if (.not. normalFieldUpdateOnly) then

     ! hy_use3dFullCTU=.true. will provide CFL <= 1; otherwise, CFL<0.5.
     if (hy_use3dFullCTU) then
        ! Set dt,dx,dy,dz factors
        dt2dxdy6=dt*dt/(6.*del(DIR_X)*del(DIR_Y))
        dt2dydz6=dt*dt/(6.*del(DIR_Y)*del(DIR_Z))
        dt2dzdx6=dt*dt/(6.*del(DIR_Z)*del(DIR_X))


        ! Now let's compute cross derivatives for CTU
        !$omp do reduction(min:minCfl) schedule(static)
        do k=k0-2,kmax+2
           do j=j0-2,jmax+2
              do i=i0-2,imax+2

#ifdef CFL_VAR
                 cellCfl = U(CFL_VAR,i,j,k)
#else
                 cellCfl = hy_cfl
#endif
                 !! save the cell center values for later use
                 Vc(HY_DENS:HY_END_VARS-kGrav) = &
                      (/U(DENS_VAR,i,j,k)&
                       ,U(VELX_VAR:VELZ_VAR,i,j,k)&
                       ,U(PRES_VAR,i,j,k)&
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
                       ,U(MAGX_VAR:MAGZ_VAR,i,j,k)&
#endif
#ifdef FLASH_UGLM_MHD
                       ,U(GLMP_VAR,i,j,k) &
#endif
                       ,U(GAMC_VAR:GAME_VAR,i,j,k) &
                       ,U(EINT_VAR,i,j,k) &
#ifdef FLASH_UHD_3T
                       ,U(EELE_VAR,i,j,k) &
                       ,U(EION_VAR,i,j,k) &
                       ,U(ERAD_VAR,i,j,k) &
#endif
                       /)

                 !! ============ X-direction ==================================================================
                 If (i .ge. i0-1 .and. i .le. imax+1) then
#ifdef FLASH_UHD_HYDRO
                 If ((j .ge. j0   .and. j .le. jmax  ) .and. (k .ge. k0 .and. k .le. kmax)) then
#endif
                 ! YZ cross dervatives for X states
                 SigmPtr => sig(:,DIR_Z,i,j-1,k)
                 SigcPtr => sig(:,DIR_Z,i,j  ,k)
                 SigpPtr => sig(:,DIR_Z,i,j+1,k)

                 call  hy_uhd_upwindTransverseFlux&
                      (dir,transOrder3D,&
                       SigmPtr,SigcPtr,SigpPtr,&
                       lambda   (1,DIR_Y,i,j,k),&
                       leftEig(1,1,DIR_Y,i,j,k),&
                       rghtEig(1,1,DIR_Y,i,j,k),&
                       HY_END_VARS,TransFluxYZ(1))

                 ! ZY cross derivatives for X states
                 SigmPtr => sig(:,DIR_Y,i,j,k-1)
                 SigcPtr => sig(:,DIR_Y,i,j,k  )
                 SigpPtr => sig(:,DIR_Y,i,j,k+1)

                 call hy_uhd_upwindTransverseFlux&
                      (dir,transOrder3D,&
                       SigmPtr,SigcPtr,SigpPtr,&
                       lambda   (1,DIR_Z,i,j,k),&
                       leftEig(1,1,DIR_Z,i,j,k),&
                       rghtEig(1,1,DIR_Z,i,j,k),&
                       HY_END_VARS,TransFluxZY(1))

#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
                 TransFluxYZ(HY_MAGX) = 0.
                 TransFluxZY(HY_MAGX) = 0.
#endif
                 scrchFaceXPtr(HY_P01_FACEXPTR_VAR:HY_P01_FACEXPTR_VAR+HY_SCRATCH_NUM-1,i,j,k)=&
                      scrchFaceXPtr(HY_P01_FACEXPTR_VAR:HY_P01_FACEXPTR_VAR+HY_SCRATCH_NUM-1,i,j,k)&
                      +(TransFluxYZ(HY_DENS:HY_END_VARS-kGrav)+TransFluxZY(HY_DENS:HY_END_VARS-kGrav))*dt2dydz6

                 scrchFaceXPtr(HY_N01_FACEXPTR_VAR:HY_N01_FACEXPTR_VAR+HY_SCRATCH_NUM-1,i,j,k)=&
                      scrchFaceXPtr(HY_N01_FACEXPTR_VAR:HY_N01_FACEXPTR_VAR+HY_SCRATCH_NUM-1,i,j,k)&
                      +(TransFluxYZ(HY_DENS:HY_END_VARS-kGrav)+TransFluxZY(HY_DENS:HY_END_VARS-kGrav))*dt2dydz6

                 !! CHECK FOR NEGATIVITY OF DENSITY AND PRESSURE IN X-DIRECTION
                 IF (scrchFaceXPtr(HY_P01_FACEXPTR_VAR,i,j,k) .le. 0. .or. &
                     scrchFaceXPtr(HY_N01_FACEXPTR_VAR,i,j,k) .le. 0. .or. &
                     scrchFaceXPtr(HY_P05_FACEXPTR_VAR,i,j,k) .le. 0. .or. &
                     scrchFaceXPtr(HY_N05_FACEXPTR_VAR,i,j,k) .le. 0. ) THEN

!!$                    scrchFaceXPtr(HY_P01_FACEXPTR_VAR:HY_P01_FACEXPTR_VAR+HY_SCRATCH_NUM-1,i,j,k)=&
!!$                         Vc(HY_DENS:HY_END_VARS-kGrav)
!!$
!!$                    scrchFaceXPtr(HY_N01_FACEXPTR_VAR:HY_N01_FACEXPTR_VAR+HY_SCRATCH_NUM-1,i,j,k)=&
!!$                         Vc(HY_DENS:HY_END_VARS-kGrav)
                    call fallbackToFirstOrder(DIR_X,&
                         scrchFaceXPtr(HY_N01_FACEXPTR_VAR:,i,j,k),&
                         scrchFaceXPtr(HY_P01_FACEXPTR_VAR:,i,j,k),&
                         Vc, &
                         hy_SpcL,hy_SpcR,U,i,j,k)
                    cellCfl = min(cellCfl, hy_cflFallbackFactor / real(NDIM))

#if (NSPECIES+NMASS_SCALARS) > 0
                 ELSE if (hy_fullSpecMsFluxHandling) then
                 ! YZ cross dervatives for X states
                    SigmPtr => hy_SpcSig(:,i,j-1,k,DIR_Z)
                    SigcPtr => hy_SpcSig(:,i,j  ,k,DIR_Z)
                    SigpPtr => hy_SpcSig(:,i,j+1,k,DIR_Z)

                    call hy_uhd_upwindTransverseFlux&
                      (dir,transOrder3D,&
                       SigmPtr,SigcPtr,SigpPtr,&
                       lambda   (1,DIR_Y,i,j,k),&
                       leftEig(1,1,DIR_Y,i,j,k),&
                       rghtEig(1,1,DIR_Y,i,j,k),&
                       HY_NSPEC,TransFluxYZ(HY_SPEC_BEG),&
                       speciesScalar=.true.)

                 ! ZY cross derivatives for X states
                    SigmPtr => hy_SpcSig(:,i,j,k-1,DIR_Y)
                    SigcPtr => hy_SpcSig(:,i,j,k  ,DIR_Y)
                    SigpPtr => hy_SpcSig(:,i,j,k+1,DIR_Y)

                    call hy_uhd_upwindTransverseFlux&
                      (dir,transOrder3D,&
                       SigmPtr,SigcPtr,SigpPtr,&
                       lambda   (1,DIR_Z,i,j,k),&
                       leftEig(1,1,DIR_Z,i,j,k),&
                       rghtEig(1,1,DIR_Z,i,j,k),&
                       HY_NSPEC,TransFluxZY(HY_SPEC_BEG),&
                       speciesScalar=.true.)

                    hy_SpcR(1:HY_NSPEC,i,j,k,DIR_X) = hy_SpcR(1:HY_NSPEC,i,j,k,DIR_X)&
                      +(TransFluxYZ(HY_SPEC_BEG:HY_SPEC_END)+TransFluxZY(HY_SPEC_BEG:HY_SPEC_END))*dt2dydz6

                    hy_SpcL(1:HY_NSPEC,i,j,k,DIR_X) = hy_SpcL(1:HY_NSPEC,i,j,k,DIR_X)&
                      +(TransFluxYZ(HY_SPEC_BEG:HY_SPEC_END)+TransFluxZY(HY_SPEC_BEG:HY_SPEC_END))*dt2dydz6
#endif
                 ENDIF
#ifdef FLASH_UHD_HYDRO
                 Endif
#endif
                 Endif

                 !! ============ y-direction ==================================================================
                 If (j .ge. j0-1 .and. j .le. jmax+1) then
#ifdef FLASH_UHD_HYDRO
                 If ((i .ge. i0   .and. i .le. imax  ) .and. (k .ge. k0 .and. k .le. kmax)) then
#endif
                 ! ZX cross derivatives for Y states
                 SigmPtr => sig(:,DIR_X,i,j,k-1)
                 SigcPtr => sig(:,DIR_X,i,j,k  )
                 SigpPtr => sig(:,DIR_X,i,j,k+1)

                 call  hy_uhd_upwindTransverseFlux&
                      (dir,transOrder3D,&
                       SigmPtr,SigcPtr,SigpPtr,&
                       lambda   (1,DIR_Z,i,j,k),&
                       leftEig(1,1,DIR_Z,i,j,k),&
                       rghtEig(1,1,DIR_Z,i,j,k),&
                       HY_END_VARS,TransFluxZX(1))

                 ! XZ cross derivatives for Y states
                 SigmPtr => sig(:,DIR_Z,i-1,j,k)
                 SigcPtr => sig(:,DIR_Z,i  ,j,k)
                 SigpPtr => sig(:,DIR_Z,i+1,j,k)

                 call  hy_uhd_upwindTransverseFlux&
                      (dir,transOrder3D,&
                       SigmPtr,SigcPtr,SigpPtr,&
                       lambda   (1,DIR_X,i,j,k),&
                       leftEig(1,1,DIR_X,i,j,k),&
                       rghtEig(1,1,DIR_X,i,j,k),&
                       HY_END_VARS,TransFluxXZ(1))

#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
                 TransFluxZX(HY_MAGY) = 0.
                 TransFluxXZ(HY_MAGY) = 0.
#endif
                 scrchFaceYPtr(HY_P01_FACEYPTR_VAR:HY_P01_FACEYPTR_VAR+HY_SCRATCH_NUM-1,i,j,k)=&
                      scrchFaceYPtr(HY_P01_FACEYPTR_VAR:HY_P01_FACEYPTR_VAR+HY_SCRATCH_NUM-1,i,j,k)&
                      + (TransFluxZX(HY_DENS:HY_END_VARS-kGrav)+TransFluxXZ(HY_DENS:HY_END_VARS-kGrav))*dt2dzdx6

                 scrchFaceYPtr(HY_N01_FACEYPTR_VAR:HY_N01_FACEYPTR_VAR+HY_SCRATCH_NUM-1,i,j,k)=&
                      scrchFaceYPtr(HY_N01_FACEYPTR_VAR:HY_N01_FACEYPTR_VAR+HY_SCRATCH_NUM-1,i,j,k)&
                      + (TransFluxZX(HY_DENS:HY_END_VARS-kGrav)+TransFluxXZ(HY_DENS:HY_END_VARS-kGrav))*dt2dzdx6

                 !! CHECK FOR NEGATIVITY OF DENSITY AND PRESSURE IN Y-DIRECTION
                 IF (scrchFaceYPtr(HY_P01_FACEYPTR_VAR,i,j,k) .le. 0. .or. &
                     scrchFaceYPtr(HY_N01_FACEYPTR_VAR,i,j,k) .le. 0. .or. &
                     scrchFaceYPtr(HY_P05_FACEYPTR_VAR,i,j,k) .le. 0. .or. &
                     scrchFaceYPtr(HY_N05_FACEYPTR_VAR,i,j,k) .le. 0. ) THEN

!!$                    scrchFaceYPtr(HY_P01_FACEYPTR_VAR:HY_P01_FACEYPTR_VAR+HY_SCRATCH_NUM-1,i,j,k)=&
!!$                         Vc(HY_DENS:HY_END_VARS-kGrav)
!!$
!!$                    scrchFaceYPtr(HY_N01_FACEYPTR_VAR:HY_N01_FACEYPTR_VAR+HY_SCRATCH_NUM-1,i,j,k)=&
!!$                         Vc(HY_DENS:HY_END_VARS-kGrav)
                    call fallbackToFirstOrder(DIR_Y,&
                         scrchFaceYPtr(HY_N01_FACEYPTR_VAR:,i,j,k),&
                         scrchFaceYPtr(HY_P01_FACEYPTR_VAR:,i,j,k),&
                         Vc, &
                         hy_SpcL,hy_SpcR,U,i,j,k)
                    cellCfl = min(cellCfl, hy_cflFallbackFactor / real(NDIM))

#if (NSPECIES+NMASS_SCALARS) > 0
                 ELSE if (hy_fullSpecMsFluxHandling) then
                 !ZX cross dervatives for X states
                    SigmPtr => hy_SpcSig(:,i,j,k-1,DIR_X)
                    SigcPtr => hy_SpcSig(:,i,j,k  ,DIR_X)
                    SigpPtr => hy_SpcSig(:,i,j,k+1,DIR_X)

                    call hy_uhd_upwindTransverseFlux&
                      (dir,transOrder3D,&
                       SigmPtr,SigcPtr,SigpPtr,&
                       lambda   (1,DIR_Z,i,j,k),&
                       leftEig(1,1,DIR_Z,i,j,k),&
                       rghtEig(1,1,DIR_Z,i,j,k),&
                       HY_NSPEC,TransFluxZX(HY_SPEC_BEG),&
                       speciesScalar=.true.)

                 ! XZ cross derivatives for Y states
                    SigmPtr => hy_SpcSig(:,i-1,j,k,DIR_Z)
                    SigcPtr => hy_SpcSig(:,i  ,j,k,DIR_Z)
                    SigpPtr => hy_SpcSig(:,i+1,j,k,DIR_Z)

                    call hy_uhd_upwindTransverseFlux&
                      (dir,transOrder3D,&
                       SigmPtr,SigcPtr,SigpPtr,&
                       lambda   (1,DIR_X,i,j,k),&
                       leftEig(1,1,DIR_X,i,j,k),&
                       rghtEig(1,1,DIR_X,i,j,k),&
                       HY_NSPEC,TransFluxXZ(HY_SPEC_BEG),&
                       speciesScalar=.true.)

                    hy_SpcR(1:HY_NSPEC,i,j,k,DIR_Y) = hy_SpcR(1:HY_NSPEC,i,j,k,DIR_Y)&
                      +(TransFluxZX(HY_SPEC_BEG:HY_SPEC_END)+TransFluxXZ(HY_SPEC_BEG:HY_SPEC_END))*dt2dydz6

                    hy_SpcL(1:HY_NSPEC,i,j,k,DIR_Y) = hy_SpcL(1:HY_NSPEC,i,j,k,DIR_Y)&
                      +(TransFluxZX(HY_SPEC_BEG:HY_SPEC_END)+TransFluxXZ(HY_SPEC_BEG:HY_SPEC_END))*dt2dydz6
#endif
                 ENDIF
#ifdef FLASH_UHD_HYDRO
                 Endif
#endif
                 Endif


                 !! ============ z-direction ==================================================================
                 If (k .ge. k0-1 .and. k .le. kmax+1) then
#ifdef FLASH_UHD_HYDRO
                 If ((i .ge. i0   .and. i .le. imax  ) .and. (j .ge. j0 .and. j .le. jmax)) then
#endif
                 ! XY cross derivatives for Z states
                 SigmPtr => sig(:,DIR_Y,i-1,j,k)
                 SigcPtr => sig(:,DIR_Y,i  ,j,k)
                 SigpPtr => sig(:,DIR_Y,i+1,j,k)

                 call  hy_uhd_upwindTransverseFlux&
                      (dir,transOrder3D,&
                       SigmPtr,SigcPtr,SigpPtr,&
                       lambda   (1,DIR_X,i,j,k),&
                       leftEig(1,1,DIR_X,i,j,k),&
                       rghtEig(1,1,DIR_X,i,j,k),&
                       HY_END_VARS,TransFluxXY(1))

                 ! YX cross derivatives for Z states
                 SigmPtr => sig(:,DIR_X,i,j-1,k)
                 SigcPtr => sig(:,DIR_X,i,j  ,k)
                 SigpPtr => sig(:,DIR_X,i,j+1,k)

                 call  hy_uhd_upwindTransverseFlux&
                      (dir,transOrder3D,&
                       SigmPtr,SigcPtr,SigpPtr,&
                       lambda   (1,DIR_Y,i,j,k),&
                       leftEig(1,1,DIR_Y,i,j,k),&
                       rghtEig(1,1,DIR_Y,i,j,k),&
                       HY_END_VARS,TransFluxYX(1))

#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
                 TransFluxXY(HY_MAGZ) = 0.
                 TransFluxYX(HY_MAGZ) = 0.
#endif

                 scrchFaceZPtr(HY_P01_FACEZPTR_VAR:HY_P01_FACEZPTR_VAR+HY_SCRATCH_NUM-1,i,j,k)=&
                      scrchFaceZPtr(HY_P01_FACEZPTR_VAR:HY_P01_FACEZPTR_VAR+HY_SCRATCH_NUM-1,i,j,k)&
                      + (TransFluxXY(HY_DENS:HY_END_VARS-kGrav)+TransFluxYX(HY_DENS:HY_END_VARS-kGrav))*dt2dxdy6

                 scrchFaceZPtr(HY_N01_FACEZPTR_VAR:HY_N01_FACEZPTR_VAR+HY_SCRATCH_NUM-1,i,j,k)=&
                      scrchFaceZPtr(HY_N01_FACEZPTR_VAR:HY_N01_FACEZPTR_VAR+HY_SCRATCH_NUM-1,i,j,k)&
                      + (TransFluxXY(HY_DENS:HY_END_VARS-kGrav)+TransFluxYX(HY_DENS:HY_END_VARS-kGrav))*dt2dxdy6


                 !! CHECK FOR NEGATIVITY OF DENSITY AND PRESSURE IN Z-DIRECTION
                 IF (scrchFaceZPtr(HY_P01_FACEZPTR_VAR,i,j,k) .le. 0. .or. &
                     scrchFaceZPtr(HY_N01_FACEZPTR_VAR,i,j,k) .le. 0. .or. &
                     scrchFaceZPtr(HY_P05_FACEZPTR_VAR,i,j,k) .le. 0. .or. &
                     scrchFaceZPtr(HY_N05_FACEZPTR_VAR,i,j,k) .le. 0. ) THEN

!!$                    scrchFaceZPtr(HY_P01_FACEZPTR_VAR:HY_P01_FACEZPTR_VAR+HY_SCRATCH_NUM-1,i,j,k)=&
!!$                         Vc(HY_DENS:HY_END_VARS-kGrav)
!!$
!!$                    scrchFaceZPtr(HY_N01_FACEZPTR_VAR:HY_N01_FACEZPTR_VAR+HY_SCRATCH_NUM-1,i,j,k)=&
!!$                         Vc(HY_DENS:HY_END_VARS-kGrav)
                    call fallbackToFirstOrder(DIR_Z,&
                         scrchFaceZPtr(HY_N01_FACEZPTR_VAR:,i,j,k),&
                         scrchFaceZPtr(HY_P01_FACEZPTR_VAR:,i,j,k),&
                         Vc, &
                         hy_SpcL,hy_SpcR,U,i,j,k)
                    cellCfl = min(cellCfl, hy_cflFallbackFactor / real(NDIM))


#if (NSPECIES+NMASS_SCALARS) > 0
                 ELSE if (hy_fullSpecMsFluxHandling) then
                 ! XY cross dervatives for X states
                    SigmPtr => hy_SpcSig(:,i-1,j,k,DIR_Y)
                    SigcPtr => hy_SpcSig(:,i,  j,k,DIR_Y)
                    SigpPtr => hy_SpcSig(:,i+1,j,k,DIR_Y)

                    call hy_uhd_upwindTransverseFlux&
                      (dir,transOrder3D,&
                       SigmPtr,SigcPtr,SigpPtr,&
                       lambda   (1,DIR_X,i,j,k),&
                       leftEig(1,1,DIR_X,i,j,k),&
                       rghtEig(1,1,DIR_X,i,j,k),&
                       HY_NSPEC,TransFluxXY(HY_SPEC_BEG),&
                       speciesScalar=.true.)

                 ! YX cross derivatives for Y states
                    SigmPtr => hy_SpcSig(:,i,j-1,k,DIR_X)
                    SigcPtr => hy_SpcSig(:,i,j,  k,DIR_X)
                    SigpPtr => hy_SpcSig(:,i,j+1,k,DIR_X)

                    call hy_uhd_upwindTransverseFlux&
                      (dir,transOrder3D,&
                       SigmPtr,SigcPtr,SigpPtr,&
                       lambda   (1,DIR_Y,i,j,k),&
                       leftEig(1,1,DIR_Y,i,j,k),&
                       rghtEig(1,1,DIR_Y,i,j,k),&
                       HY_NSPEC,TransFluxYX(HY_SPEC_BEG),&
                       speciesScalar=.true.)

                    hy_SpcR(1:HY_NSPEC,i,j,k,DIR_Z) = hy_SpcR(1:HY_NSPEC,i,j,k,DIR_Z)&
                      +(TransFluxXY(HY_SPEC_BEG:HY_SPEC_END)+TransFluxYX(HY_SPEC_BEG:HY_SPEC_END))*dt2dydz6

                    hy_SpcL(1:HY_NSPEC,i,j,k,DIR_Z) = hy_SpcL(1:HY_NSPEC,i,j,k,DIR_Z)&
                      +(TransFluxXY(HY_SPEC_BEG:HY_SPEC_END)+TransFluxYX(HY_SPEC_BEG:HY_SPEC_END))*dt2dydz6
#endif
                 ENDIF
#ifdef FLASH_UHD_HYDRO
                 Endif
#endif
                 Endif

                 if (hy_fallbackLowerCFL) then
#ifdef CFL_VAR
                    U(CFL_VAR,i,j,k) = cellCfl
#else
                    minCfl = min(minCfl,cellCfl)
#endif
                 end if

              enddo ! i-loop
           enddo ! j-loop
        enddo ! k-loop
        !$omp end do

#ifndef CFL_VAR
        if (hy_fallbackLowerCFL) then
           !$omp critical (hy_crit_update_hy_cfl)
           hy_cfl = min(hy_cfl,minCfl)
           !$omp end critical (hy_crit_update_hy_cfl)
        end if
#endif

     end if !end of if (hy_use3dFullCTU) then

  endif ! End of if (.not. normalFieldUpdateOnly) then

#endif /* end of #if NDIM == 3 */

  !$omp end parallel




!#ifdef GLMGLM
#ifdef FLASH_UGLM_MHD
  do k=k0-2-k3+kHydro*k3,kmax+2+k3-kHydro*k3
     do j=j0-2-k3+kHydro*k2,jmax+2+k3-kHydro*k2
        do i=i0-2-k3+kHydro,imax+2+k3-kHydro
           ! Extra stencil is needed for 3D to correctly calculate transverse fluxes 
           !(i.e., cross derivatives in x,y, & z)

           ! (1) costruct Godunov flux for GLM-Psi (GMLP_VAR)
           BxStar(i,j,k) = ( scrchFaceXPtr(HY_P06_FACEXPTR_VAR,i-1,j,  k  )&
                            +scrchFaceXPtr(HY_N06_FACEXPTR_VAR,i,  j,  k  ))*0.5

           BxStar(i,j,k) = BxStar(i,j,k) - 0.5/hy_C_hyp*&
                           ( scrchFaceXPtr(HY_N09_FACEXPTR_VAR,i,  j,  k  ) &
                            -scrchFaceXPtr(HY_P09_FACEXPTR_VAR,i-1,j,  k  ))
#if NDIM > 1
           ByStar(i,j,k) = ( scrchFaceYPtr(HY_P07_FACEYPTR_VAR,i,  j-1,k  )&
                            +scrchFaceYPtr(HY_N07_FACEYPTR_VAR,i,  j,  k  ))*0.5

           ByStar(i,j,k) = ByStar(i,j,k) - 0.5/hy_C_hyp*&
                           ( scrchFaceYPtr(HY_N09_FACEYPTR_VAR,i,  j,  k  ) &
                            -scrchFaceYPtr(HY_P09_FACEYPTR_VAR,i,  j-1,k  ))
#if NDIM == 3
           BzStar(i,j,k) = ( scrchFaceZPtr(HY_P08_FACEZPTR_VAR,i,  j,  k-1)&
                            +scrchFaceZPtr(HY_N08_FACEZPTR_VAR,i,  j,  k  ))*0.5

           BzStar(i,j,k) = BzStar(i,j,k) - 0.5/hy_C_hyp*&
                           ( scrchFaceZPtr(HY_N09_FACEZPTR_VAR,i,  j,  k  ) &
                            -scrchFaceZPtr(HY_P09_FACEZPTR_VAR,i,  j,  k-1))
#endif
#endif


           ! (2) costruct Godunov flux for normal magnetic fields
           GLMxStar(i,j,k) = ( scrchFaceXPtr(HY_P09_FACEXPTR_VAR,i-1,j,  k  )&
                              +scrchFaceXPtr(HY_N09_FACEXPTR_VAR,i,  j,  k  ))*0.5

           GLMxStar(i,j,k) = GLMxStar(i,j,k) - 0.5*hy_C_hyp*&
                             ( scrchFaceXPtr(HY_N06_FACEXPTR_VAR,i,  j,  k  ) &
                              -scrchFaceXPtr(HY_P06_FACEXPTR_VAR,i-1,j,  k  ))
#if NDIM > 1
           GLMyStar(i,j,k) = ( scrchFaceYPtr(HY_P09_FACEYPTR_VAR,i,  j-1,k  )&
                              +scrchFaceYPtr(HY_N09_FACEYPTR_VAR,i,  j,  k  ))*0.5

           GLMyStar(i,j,k) = GLMyStar(i,j,k) - 0.5*hy_C_hyp*&
                             ( scrchFaceYPtr(HY_N07_FACEYPTR_VAR,i,  j,  k  ) &
                              -scrchFaceYPtr(HY_P07_FACEYPTR_VAR,i,  j-1,k  ))
#if NDIM == 3
           GLMzStar(i,j,k) = ( scrchFaceZPtr(HY_P09_FACEZPTR_VAR,i,  j,  k-1)&
                              +scrchFaceZPtr(HY_N09_FACEZPTR_VAR,i,  j,  k  ))*0.5

           GLMzStar(i,j,k) = GLMzStar(i,j,k) - 0.5*hy_C_hyp*&
                             ( scrchFaceZPtr(HY_N08_FACEZPTR_VAR,i,  j,  k  ) &
                              -scrchFaceZPtr(HY_P08_FACEZPTR_VAR,i,  j,  k-1))
#endif
#endif

        enddo
     enddo
  enddo


  do k=k0-2-k3+kHydro*k3,kmax+2+k3-kHydro*k3
     do j=j0-2-k3+kHydro*k2,jmax+2+k3-kHydro*k2
        do i=i0-2-k3+kHydro,imax+2+k3-kHydro
           ! GLM psi x-Riemann state
           scrchFaceXPtr(HY_P09_FACEXPTR_VAR,i,j,k) = GLMxSTar(i+1,j,  k  )
           scrchFaceXPtr(HY_N09_FACEXPTR_VAR,i,j,k) = GLMxStar(i,  j,  k  )

           ! magx
           scrchFaceXPtr(HY_P06_FACEXPTR_VAR,i,j,k) = BxStar  (i+1,j,  k  )
           scrchFaceXPtr(HY_N06_FACEXPTR_VAR,i,j,k) = BxStar  (i,  j,  k  )

#if NDIM > 1
           ! GLM psi y-Riemann state
           scrchFaceYPtr(HY_P09_FACEYPTR_VAR,i,j,k) = GLMySTar(i,  j+1,k  )
           scrchFaceYPtr(HY_N09_FACEYPTR_VAR,i,j,k) = GLMyStar(i,  j,  k  )

           scrchFaceYPtr(HY_P07_FACEYPTR_VAR,i,j,k) = ByStar  (i,  j+1,k  )
           scrchFaceYPtr(HY_N07_FACEYPTR_VAR,i,j,k) = ByStar  (i,  j,  k  )

#if NDIM == 3
           ! GLM psi z-Riemann state
           scrchFaceZPtr(HY_P09_FACEZPTR_VAR,i,j,k) = GLMzSTar(i,  j,  k+1)
           scrchFaceZPtr(HY_N09_FACEZPTR_VAR,i,j,k) = GLMzStar(i,  j,  k  )

           scrchFaceZPtr(HY_P08_FACEZPTR_VAR,i,j,k) = BzStar  (i,  j,  k+1)
           scrchFaceZPtr(HY_N08_FACEZPTR_VAR,i,j,k) = BzStar  (i,  j,  k  )
#endif
#endif

        enddo
     enddo
  enddo

#endif /* ifdef FLASH_UGLM_MHD */
!#endif


  !! Release pointers
  call Grid_releaseBlkPtr(blockID,U,CENTER)

  ! MHD only-------------------------------------------------------------------------------
#if defined(FLASH_USM_MHD) && NFACE_VARS > 0 && NDIM > 1
  if (hy_order > 1) then
     call Grid_releaseBlkPtr(blockID,Bx,FACEX)
     call Grid_releaseBlkPtr(blockID,By,FACEY)
     if (NDIM == 3) call Grid_releaseBlkPtr(blockID,Bz,FACEZ)
  endif ! if (hy_order > 1) then
#endif /* endif of if defined(FLASH_USM_MHD) && NFACE_VARS > 0 && NDIM > 1 */
  ! MHD only-------------------------------------------------------------------------------


  !! Deallocate arrays
  !deallocate(DivU)
  if (.NOT. normalFieldUpdateOnly) then
     deallocate(sig)
     deallocate(lambda)
     deallocate(leftEig)
     deallocate(rghtEig)
  end if

contains
#include "FortranLangFeatures.fh"
  subroutine fallbackToFirstOrder(iDir,Wleft,Wright,Vc,spcL,spcR,U,i,j,k)
    integer,intent(IN)                        :: iDir
    real,   intent(OUT), dimension(:)         :: Wleft,Wright
    real,   intent(IN),  dimension(:)         :: Vc
    real,   POINTER_INTENT_IN, dimension(:,:,:,:,:),OPTIONAL :: spcL,spcR
    real,   intent(IN),  dimension(:,:,:,:)  ,OPTIONAL :: U
    integer,intent(IN),                       OPTIONAL :: i,j,k

    Wleft (HY_DENS:HY_END_VARS-kGrav) = Vc(HY_DENS:HY_END_VARS-kGrav)
    Wright(HY_DENS:HY_END_VARS-kGrav) = Vc(HY_DENS:HY_END_VARS-kGrav)
    if (hy_fullSpecMsFluxHandling .AND. hy_numXN > 0 &
         .AND. present(spcR)) then
       spcL(:,i,j,k,iDir) = &
            U(SPECIES_BEGIN:MASS_SCALARS_END,i,j,k)
       spcR(:,i,j,k,iDir) = &
            U(SPECIES_BEGIN:MASS_SCALARS_END,i,j,k)
    end if

  end subroutine fallbackToFirstOrder
End Subroutine hy_uhd_getRiemannState
