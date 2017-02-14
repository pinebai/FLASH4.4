!!****if* source/physics/Hydro/HydroMain/unsplit/hy_uhd_DataReconstructNormalDir_PPM
!!
!! NAME
!!
!!  hy_uhd_DataReconstructNormalDir_PPM
!!
!! SYNOPSIS
!!
!!  call hy_uhd_DataReconstructNormalDir_PPM(integer(IN) :: dir,
!!                                    real(IN)    :: dt,
!!                                    real(IN)    :: delta,
!!                                    pointer(IN) :: Data1D(:,:),
!!                                    pointer(IN) :: DataGrav1D(:),
!!                                    real(IN)    :: FlatCoeff,
!!                                    logical(IN) :: TransUpdateOnly,
!!                                    real(OUT)   :: lambda0(:),
!!                                    real(OUT)   :: leig0(:,:),
!!                                    real(OUT)   :: reig0(:,:),
!!                                    real(OUT)   :: Wp(:),
!!                                    real(OUT)   :: Wm(:),
!!                                    real(OUT)   :: sig(:),
!!                                    real(IN),optional :: dnBnormal,
!!                                    real(IN),optional :: dnGLMnormal,
!!                                    real(IN),optional :: aBnormal(:),
!!                                    real(IN),optional :: aGLMnormal(:),
!!                                    real(OUT),optional :: Sr(:),
!!                                    real(OUT),optional :: Sl(:),
!!                                    real(OUT),optional :: SpcSig(:),
!!
!! ARGUMENTS
!!
!!  dir         - normal direction along which the reconstuction is performed
!!  dt          - timestep
!!  delta       - deltas in each {x,y,z} direction
!!  Data1D      - pointer array holding neighboring stencil hydro/MHD data for reconstruction
!!  DataGrav1D  - pointer array holding neighboring stencil gravity data for reconstruction
!!  FlatCoeff   - flattening parameters, primarily for PPM
!!  TransUpdateOnly - a switch for a selective transverse flux update in the normal direction
!!  lambda      - eigenvalue
!!  leig        - left eigenvector
!!  reig        - right eigenvector
!!  Wp,Wm       - left(minus) and right(plus) Riemann states
!!  sig         - transverse flux term
!!  dnBnormal   - MHD multidimensional term for constrained-transport (i.e., USM) MHD scheme
!!  dnGLMnormal - additional term for GLM-MHD
!!  aBnormal    - MHD multidimensional term for constrained-transport (i.e., USM) MHD scheme
!!  aGLMnormal  - additional term for GLM-MHD
!!  Sr          - right Riemann state of species and mass scalars 
!!  Sl          - left  Riemann state of species and mass scalars 
!!  SpcSig      - transverse flux term for spieces and mass scalars
!!
!!
!! DESCRIPTION
!!
!!  This subroutine provides a third-order spatially accurate piecewise parabolic method (PPM) for
!!  reconstruction.
!!
!! REFERENCES
!!
!!  * Colella and Woodward, 54, 174 (1984), JCP
!!  * Toro, Riemann Solvers and Numerical Methods for Fluid Dynamics, Springer, 1997
!!  * Lee, D., Ph.D. Dissertation, Univ. of MD, 2006
!!  * Lee, D. and Deane, A., "An Unsplit Staggered Mesh Scheme for Multidimensional
!!                            Magnetohydrodynamics", 228 (2009), 952-975, JCP
!!  * Stone, Gardiner, Teuben, Hawley, Simon, "Athena: A new code for astrophysical MHD"
!!    arXiv:0804.0402v1 [astro-ph] 2 Apr 2008
!!  * Colella, 87, 171-200 (1990), JCP
!!
!!***

Subroutine hy_uhd_DataReconstructNormalDir_PPM&
     (dir,dt,delta,Data1D,DataGrav1D,&
      FlatCoeff,TransUpdateOnly, &
      lambda0,leig0,reig0,&
      Wp,Wm,sig,&
      dnBnormal,aBnormal,&     !These are optional
      dnGLMnormal,aGLMnormal,& !These are optional
      Sr,Sl,SpcSig)            !These are optional

  use Hydro_data,           ONLY : hy_charLimiting,   &
                                   hy_eswitch,        &
                                   hy_RiemannSolver,  &
                                   hy_tiny,           &
                                   hy_entropy,        &
                                   hy_flattening,     &
                                   hy_transOrder,     &
                                   hy_ContactSteepening,&
                                   hy_upwindTVD,      &
                                   hy_3Torder,        &
                                   hy_useAuxEintEqn,  &
                                   hy_fullSpecMsFluxHandling

  use hy_uhd_interface,     ONLY : hy_uhd_TVDslope,       &
                                   hy_uhd_TVDslopeUpwind, & 
                                   hy_uhd_upwindTransverseFlux,&
                                   hy_uhd_eigenParameters, &
                                   hy_uhd_eigenValue,      &
                                   hy_uhd_eigenVector

  use hy_uhd_slopeLimiters, ONLY : minmod
  use Timers_interface,     ONLY : Timers_start, Timers_stop

  implicit none

#include "constants.h"
#include "Flash.h"
#include "UHD.h"

  !!-----Arguments---------------------------------------------------------
  integer,intent(IN) :: dir
  real,   intent(IN) :: dt,delta
  real, pointer, dimension(:,:) :: Data1D
  real, pointer, dimension(:)   :: DataGrav1D
  real,    intent(IN) :: FlatCoeff
  logical, intent(IN) :: TransUpdateOnly
  real, dimension(HY_WAVENUM),intent(OUT) :: lambda0
  real, dimension(HY_VARINUM,HY_WAVENUM),intent(OUT) :: leig0
  real, dimension(HY_VARINUM,HY_WAVENUM),intent(OUT) :: reig0
  real,intent(OUT),dimension(HY_VARINUMMAX) :: Wp,Wm
  real,intent(OUT),dimension(HY_VARINUMMAX) :: sig
  !optional arguments for MHD-----
  real,intent(IN), optional :: dnBnormal,dnGLMnormal
  real,intent(IN), dimension(HY_VARINUM), optional :: aBnormal,aGLMnormal
  !-------------------------------
  real,intent(OUT),dimension(HY_NSPEC),   optional :: Sr,Sl
  real,intent(OUT),dimension(:),          optional :: SpcSig
  !!------------------------------------------------------------------------

  integer :: n, nVar, iBeg, hyEndVar,hyEndPrimVar
  real    :: dtn,hdtn,qdtn,factor
  real    :: constA,constB,constC,constD,lambdaMax,lambdaMin
  real    :: temp1, temp2, temp3, eta_steep, del2rhoR, del2rhoL, Flattening
  real, dimension(HY_VARINUMMAX) :: sigL,sigR
  real, dimension(HY_VARINUMMAX) :: vecL,vecR,delW,W6,delbar0,delbarP,delbarN
  real, dimension(HY_VARINUMMAX),target :: Vmmm,Vmm,Vm,Vc,Vp,Vpp,Vppp
  real, PARAMETER :: eta1=20.E0, eta2=0.05E0,epsln=0.01E0,K0=0.1E0

  ! EIG SYSTEM
  logical :: cons=.false.
  real    :: cf,uN,cs,ca,as,af,hyp
  real, pointer, dimension(:) :: vm2Ptr,vm1Ptr,vc0Ptr,vp1Ptr,vp2Ptr
  real, dimension(MDIM) :: beta
  real, dimension(HY_WAVENUM)   :: lambdaP,lambdaN,lambdaPP,lambdaNN
  real, dimension(HY_VARINUM,HY_WAVENUM) :: leigP,leigN,reigP,reigN


#if (NSPECIES+NMASS_SCALARS) > 0
  ! Species and mass scalars
  integer :: isph, ispu !isph for hydro scope index, ispu for unk index
  real, dimension(HY_NSPEC), target :: Sc,Sp,Sm,Spp,Smm
  real :: delbarSp0,delbarSpP,delbarSpN,Sp6
#endif


  ! Debugging mode compilation complains when these are not initialized.
  if (present(Sr)) Sr = 0.
  if (present(Sl)) Sl = 0.
  if (present(SpcSig)) SpcSig = 0.


  ! Set index range depending on hydro or MHD
  ! default for for hydro
  hyEndVar  = HY_ENER
#ifdef FLASH_USM_MHD /* for USM-MHD */
  hyEndVar  = HY_MAGZ
#elif defined(FLASH_UGLM_MHD)
  hyEndVar  = HY_MAGZ
#endif

  hyEndPrimVar = hyEndVar
  if (hy_useAuxEintEqn) hyEndVar = HY_EINT
#ifdef FLASH_UHD_3T
  hyEndVar = HY_END_VARS
#endif
#ifdef GRAVITY
  hyEndVar = HY_GRAV
#endif

  iBeg = 0
  If (TransUpdateOnly) Then
     iBeg = iBeg-1
  Endif
  !! ----------------------------------------------------------------------
  !! 1D array under consideration
  !! ----------------------------------------------------------------------

  !! [A] Usual hydro variables
  IF (.not. TransUpdateOnly) THEN
#if NGUARD > 4
     iBeg = 1

  !Always initialize with zero before storing arrays!
     Vmmm=0; Vppp=0.
     Vmmm(HY_DENS)         = Data1D(DENS_VAR,           iBeg)
     Vmmm(HY_VELX:HY_VELZ) = Data1D(VELX_VAR:VELZ_VAR,  iBeg)
     Vmmm(HY_PRES)         = Data1D(PRES_VAR,           iBeg)
     Vmmm(HY_GAMC:HY_GAME) = Data1D(GAMC_VAR:GAME_VAR,  iBeg)
     Vmmm(HY_EINT)         = Data1D(EINT_VAR,           iBeg)

     Vppp(HY_DENS)         = Data1D(DENS_VAR,         6+iBeg)
     Vppp(HY_VELX:HY_VELZ) = Data1D(VELX_VAR:VELZ_VAR,6+iBeg)
     Vppp(HY_PRES)         = Data1D(PRES_VAR,         6+iBeg)
     Vppp(HY_GAMC:HY_GAME) = Data1D(GAMC_VAR:GAME_VAR,6+iBeg)
     Vppp(HY_EINT)         = Data1D(EINT_VAR,         6+iBeg)
#endif

     Vmm=0.;Vpp=0.
     Vmm(HY_DENS)          = Data1D(DENS_VAR,         1+iBeg)
     Vmm(HY_VELX:HY_VELZ)  = Data1D(VELX_VAR:VELZ_VAR,1+iBeg)
     Vmm(HY_PRES)          = Data1D(PRES_VAR,         1+iBeg)
     Vmm(HY_GAMC:HY_GAME)  = Data1D(GAMC_VAR:GAME_VAR,1+iBeg)
     Vmm(HY_EINT)          = Data1D(EINT_VAR,         1+iBeg)

     Vpp(HY_DENS)          = Data1D(DENS_VAR,         5+iBeg)
     Vpp(HY_VELX:HY_VELZ)  = Data1D(VELX_VAR:VELZ_VAR,5+iBeg)
     Vpp(HY_PRES)          = Data1D(PRES_VAR,         5+iBeg)
     Vpp(HY_GAMC:HY_GAME)  = Data1D(GAMC_VAR:GAME_VAR,5+iBeg)
     Vpp(HY_EINT)          = Data1D(EINT_VAR,         5+iBeg)
  ENDIF

  Vm=0.;Vc=0.;Vp=0.
  Vm (HY_DENS)          = Data1D(DENS_VAR,         2+iBeg)
  Vm (HY_VELX:HY_VELZ)  = Data1D(VELX_VAR:VELZ_VAR,2+iBeg)
  Vm (HY_PRES)          = Data1D(PRES_VAR,         2+iBeg)
  Vm (HY_GAMC:HY_GAME)  = Data1D(GAMC_VAR:GAME_VAR,2+iBeg)
  Vm (HY_EINT)          = Data1D(EINT_VAR,         2+iBeg)

  Vc (HY_DENS)          = Data1D(DENS_VAR,         3+iBeg)
  Vc (HY_VELX:HY_VELZ)  = Data1D(VELX_VAR:VELZ_VAR,3+iBeg)
  Vc (HY_PRES)          = Data1D(PRES_VAR,         3+iBeg)
  Vc (HY_GAMC:HY_GAME)  = Data1D(GAMC_VAR:GAME_VAR,3+iBeg)
  Vc (HY_EINT)          = Data1D(EINT_VAR,         3+iBeg)

  Vp (HY_DENS)          = Data1D(DENS_VAR,         4+iBeg)
  Vp (HY_VELX:HY_VELZ)  = Data1D(VELX_VAR:VELZ_VAR,4+iBeg)
  Vp (HY_PRES)          = Data1D(PRES_VAR,         4+iBeg)
  Vp (HY_GAMC:HY_GAME)  = Data1D(GAMC_VAR:GAME_VAR,4+iBeg)
  Vp (HY_EINT)          = Data1D(EINT_VAR,         4+iBeg)



  !! [B] magnetic fields for MHD
#ifdef FLASH_USM_MHD
  IF (.not. TransUpdateOnly) THEN
     Vmm (HY_MAGX:HY_MAGZ) = Data1D(MAGX_VAR:MAGZ_VAR,1+iBeg)
     Vpp (HY_MAGX:HY_MAGZ) = Data1D(MAGX_VAR:MAGZ_VAR,5+iBeg)
#if NGUARD > 4
     Vmmm(HY_MAGX:HY_MAGZ) = Data1D(MAGX_VAR:MAGZ_VAR,  iBeg)
     Vppp(HY_MAGX:HY_MAGZ) = Data1D(MAGX_VAR:MAGZ_VAR,6+iBeg)
#endif
  ENDIF
  Vm  (HY_MAGX:HY_MAGZ) = Data1D(MAGX_VAR:MAGZ_VAR,2+iBeg)
  Vc  (HY_MAGX:HY_MAGZ) = Data1D(MAGX_VAR:MAGZ_VAR,3+iBeg)
  Vp  (HY_MAGX:HY_MAGZ) = Data1D(MAGX_VAR:MAGZ_VAR,4+iBeg)
#endif

  !! [C] 3T variables
#ifdef FLASH_UHD_3T
  IF (.not. TransUpdateOnly) THEN
     Vmm(HY_EELE:HY_ERAD)  = (/Data1D(EELE_VAR,1+iBeg),Data1D(EION_VAR,1+iBeg),Data1D(ERAD_VAR,1+iBeg)/)
     Vpp(HY_EELE:HY_ERAD)  = (/Data1D(EELE_VAR,5+iBeg),Data1D(EION_VAR,5+iBeg),Data1D(ERAD_VAR,5+iBeg)/)
#if NGUARD > 4
     Vmmm(HY_EELE:HY_ERAD) = (/Data1D(EELE_VAR,  iBeg),Data1D(EION_VAR,  iBeg),Data1D(ERAD_VAR,  iBeg)/)
     Vppp(HY_EELE:HY_ERAD) = (/Data1D(EELE_VAR,6+iBeg),Data1D(EION_VAR,6+iBeg),Data1D(ERAD_VAR,6+iBeg)/)
#endif
  Endif
  Vm (HY_EELE:HY_ERAD)  = (/Data1D(EELE_VAR,2+iBeg),Data1D(EION_VAR,2+iBeg),Data1D(ERAD_VAR,2+iBeg)/)
  Vc (HY_EELE:HY_ERAD)  = (/Data1D(EELE_VAR,3+iBeg),Data1D(EION_VAR,3+iBeg),Data1D(ERAD_VAR,3+iBeg)/)
  Vp (HY_EELE:HY_ERAD)  = (/Data1D(EELE_VAR,4+iBeg),Data1D(EION_VAR,4+iBeg),Data1D(ERAD_VAR,4+iBeg)/)
#endif

  !! [D] Gravity component
#ifdef GRAVITY
  IF (.not. TransUpdateOnly) THEN
     Vmm(HY_GRAV)  = DataGrav1D(1+iBeg)
     Vpp(HY_GRAV)  = DataGrav1D(5+iBeg)
#if NGUARD > 4
     Vmmm(HY_GRAV) = DataGrav1D(  iBeg)
     Vppp(HY_GRAV) = DataGrav1D(6+iBeg)
#endif
  Endif
  Vm (HY_GRAV)  = DataGrav1D(2+iBeg)
  Vc (HY_GRAV)  = DataGrav1D(3+iBeg)
  Vp (HY_GRAV)  = DataGrav1D(4+iBeg)
#endif

  !! [E] Species
#if (NSPECIES+NMASS_SCALARS) > 0
  if (hy_fullSpecMsFluxHandling) then
     do ispu = SPECIES_BEGIN, MASS_SCALARS_END !SPECIES_END
        isph= ispu-NPROP_VARS
        IF (.not. TransUpdateOnly) THEN
           Smm(isph)  = Data1D(ispu,1+iBeg)
           Spp(isph)  = Data1D(ispu,5+iBeg)
        Endif
        Sm (isph)  = Data1D(ispu,2+iBeg)
        Sc (isph)  = Data1D(ispu,3+iBeg)
        Sp (isph)  = Data1D(ispu,4+iBeg)
     enddo
  endif
#endif


  !!**********************!!
  !! BEGIN RECONSTRUCTION !!
  !!**********************!!

  !! half dt & delta
  dtn=dt/delta
  hdtn=0.5*dtn
  qdtn=0.25*dtn

  !! initialize arrays with zero
  Wp=0.;      Wm=0.
!!$ delbar0=0.; delbarN=0.; delbarP=0.

  call hy_uhd_eigenParameters&
       (Vc(HY_DENS:HY_GAME),dir,uN,cf&
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
       ,C_alfn=ca,C_slow=cs,A_f=af,A_s=as,B_beta=beta&
#ifdef FLASH_UGLM_MHD
       ,C_hyp=hyp&
#endif
#endif
       )
  call hy_uhd_eigenValue&
       (lambda0,uN,cf&
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
       ,C_alfn=ca,C_slow=cs&
#ifdef FLASH_UGLM_MHD
       ,C_hyp=hyp&
#endif
#endif
       )
  call hy_uhd_eigenVector&
       (leig0,reig0,Vc(HY_DENS:HY_GAME),dir,.false.,cf&
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
       ,C_alfn=ca,C_slow=cs,A_f=af,A_s=as,B_beta=beta&
#ifdef FLASH_UGLM_MHD
       ,C_hyp=hyp&
#endif
#endif
       )

  !! -------------------------------------------------------------------------------------!
  !! [1] First we compute transverse flux if we are interested in running ----------------!
  !!     multi-dimensional cases ---------------------------------------------------------!
  !! -------------------------------------------------------------------------------------!
  if (NDIM > 1) then

     vm1Ptr => Vm (:)
     vc0Ptr => Vc (:)
     vp1Ptr => Vp (:)

     call hy_uhd_upwindTransverseFlux&
          (dir,hy_transOrder,vm1Ptr,vc0Ptr,vp1Ptr,lambda0,leig0,reig0,HY_END_VARS,sig)


#if (NSPECIES+NMASS_SCALARS) > 0
     if (hy_fullSpecMsFluxHandling .AND. present(SpcSig)) then
        vm1Ptr => Sm (:)
        vc0Ptr => Sc (:)
        vp1Ptr => Sp (:)

        call hy_uhd_upwindTransverseFlux&
          (dir,hy_transOrder,vm1Ptr,vc0Ptr,vp1Ptr,lambda0,leig0,reig0,HY_NSPEC,SpcSig,speciesScalar=.true.)
     endif ! (hy_fullSpecMsFluxHandling)
#endif
  endif ! NDIM > 1



  !! -------------------------------------------------------------------------------------!
  !! [2] Apply TVD slope limiter for normal gradients ------------------------------------!
  !! -------------------------------------------------------------------------------------!
  IF (.not. TransUpdateOnly) THEN

  call hy_uhd_eigenParameters&
       (Vp(HY_DENS:HY_GAME),dir,uN,cf&
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
       ,C_alfn=ca,C_slow=cs,A_f=af,A_s=as,B_beta=beta&
#ifdef FLASH_UGLM_MHD
       ,C_hyp=hyp&
#endif
#endif
       )
  call hy_uhd_eigenValue&
       (lambdaP,uN,cf&
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
       ,C_alfn=ca,C_slow=cs&
#ifdef FLASH_UGLM_MHD
       ,C_hyp=hyp&
#endif
#endif
       )
  call hy_uhd_eigenVector&
       (leigP,reigP,Vp(HY_DENS:HY_GAME),dir,.false.,cf&
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
       ,C_alfn=ca,C_slow=cs,A_f=af,A_s=as,B_beta=beta&
#ifdef FLASH_UGLM_MHD
       ,C_hyp=hyp&
#endif
#endif
       )


  call hy_uhd_eigenParameters&
       (Vm(HY_DENS:HY_GAME),dir,uN,cf&
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
       ,C_alfn=ca,C_slow=cs,A_f=af,A_s=as,B_beta=beta&
#ifdef FLASH_UGLM_MHD
       ,C_hyp=hyp&
#endif
#endif
       )
  call hy_uhd_eigenValue&
       (lambdaN,uN,cf&
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
       ,C_alfn=ca,C_slow=cs&
#ifdef FLASH_UGLM_MHD
       ,C_hyp=hyp&
#endif
#endif
       )
  call hy_uhd_eigenVector&
       (leigN,reigN,Vm(HY_DENS:HY_GAME),dir,.false.,cf&
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
       ,C_alfn=ca,C_slow=cs,A_f=af,A_s=as,B_beta=beta&
#ifdef FLASH_UGLM_MHD
       ,C_hyp=hyp&
#endif
#endif
       )


     if (.not. hy_upwindTVD) then        
        call hy_uhd_TVDslope(dir,Vmm,Vm, Vc, lambdaN,leigN,delbarN)
        call hy_uhd_TVDslope(dir,Vm, Vc, Vp, lambda0,leig0,delbar0)
        call hy_uhd_TVDslope(dir,Vc, Vp, Vpp,lambdaP,leigP,delbarP)
     else


        !! lambdaNN
        call hy_uhd_eigenParameters&
             (Vmm(HY_DENS:HY_GAME),dir,uN,cf&
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
             ,C_alfn=ca,C_slow=cs,A_f=af,A_s=as,B_beta=beta&
#ifdef FLASH_UGLM_MHD
             ,C_hyp=hyp&
#endif
#endif
             )
        call hy_uhd_eigenValue&
             (lambdaNN,uN,cf&
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
             ,C_alfn=ca,C_slow=cs&
#ifdef FLASH_UGLM_MHD
             ,C_hyp=hyp&
#endif
#endif
             )

        !! lambdaPP
        call hy_uhd_eigenParameters&
             (Vpp(HY_DENS:HY_GAME),dir,uN,cf&
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
             ,C_alfn=ca,C_slow=cs,A_f=af,A_s=as,B_beta=beta&
#ifdef FLASH_UGLM_MHD
             ,C_hyp=hyp&
#endif
#endif
             )
        call hy_uhd_eigenValue&
             (lambdaPP,uN,cf&
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
             ,C_alfn=ca,C_slow=cs&
#ifdef FLASH_UGLM_MHD
             ,C_hyp=hyp&
#endif
#endif
             )


        !! Upwinded TVD for PPM ---"
        if ((lambdaN(HY_FASTRGHT) > 0. .and. lambdaP(HY_FASTLEFT) < 0.) .or. &
            (lambdaN(HY_FASTRGHT) > 0. .and. lambda0(HY_FASTLEFT) < 0.) .or. &
            (lambda0(HY_FASTRGHT) > 0. .and. lambdaP(HY_FASTLEFT) < 0.)) then
           call hy_uhd_TVDslopeUpwind&
                (dir,Vmm,Vm,Vc,Vp,Vpp,lambdaN,lambda0,lambdaP,leig0,delbar0)
        else
           call hy_uhd_TVDslope(dir,Vm,Vc,Vp,lambda0,leig0,delbar0)
        endif

        !! Upwinded TVD for PPM ---
        if ((lambda0(HY_FASTRGHT) > 0. .and. lambdaPP(HY_FASTLEFT) < 0.) .or. &
            (lambda0(HY_FASTRGHT) > 0. .and. lambdaP (HY_FASTLEFT) < 0.) .or. &
            (lambdaP(HY_FASTRGHT) > 0. .and. lambdaPP(HY_FASTLEFT) < 0.)) then
           call hy_uhd_TVDslopeUpwind&
                (dir,Vm,Vc,Vp,Vpp,Vppp,lambda0,lambdaP,lambdaPP,leigP,delbarP)
        else
           call hy_uhd_TVDslope(dir,Vc,Vp,Vpp,lambdaP,leigP,delbarP)
        endif

        !! Upwinded TVD for PPM ---
        if ((lambdaNN(HY_FASTRGHT) > 0. .and. lambda0(HY_FASTLEFT) < 0.) .or. &
            (lambdaNN(HY_FASTRGHT) > 0. .and. lambdaN(HY_FASTLEFT) < 0.) .or. &
            (lambdaN (HY_FASTRGHT) > 0. .and. lambda0(HY_FASTLEFT) < 0.)) then
           call hy_uhd_TVDslopeUpwind&
                (dir,Vmmm,Vmm,Vm,Vc,Vp,lambdaNN,lambdaN,lambda0,leigN,delbarN)
        else
           call hy_uhd_TVDslope(dir,Vmm,Vm,Vc,lambdaN,leigN,delbarN)
        endif

     endif !end of if (.not. hy_upwindTVD) then


     !! First initialize flattening coefficients
     if (hy_flattening) then
        Flattening = FlatCoeff
     else
        Flattening = 0.
     endif

     ! delbar contains characteristic variables that are limited.
     ! Project delbar to primitive variables now.
     if (hy_charLimiting) then
        ! delbar at i -----------------------------------------------------------------------
        delbar0(HY_DENS:hyEndPrimVar)  = &
             reig0(HY_DENS:hyEndPrimVar,HY_FASTLEFT)*delbar0(HY_FASTLEFT)+&
#ifdef FLASH_USM_MHD
             reig0(HY_DENS:hyEndPrimVar,HY_ALFNLEFT)*delbar0(HY_ALFNLEFT)+&
#endif
             reig0(HY_DENS:hyEndPrimVar,HY_SLOWLEFT)*delbar0(HY_SLOWLEFT)+&
             reig0(HY_DENS:hyEndPrimVar,HY_ENTROPY )*delbar0(HY_ENTROPY )+&
             reig0(HY_DENS:hyEndPrimVar,HY_SLOWRGHT)*delbar0(HY_SLOWRGHT)+&
#ifdef FLASH_USM_MHD
             reig0(HY_DENS:hyEndPrimVar,HY_ALFNRGHT)*delbar0(HY_ALFNRGHT)+&
#endif
             reig0(HY_DENS:hyEndPrimVar,HY_FASTRGHT)*delbar0(HY_FASTRGHT)

        ! delbar at i+1 -------------------------------------------------------------------------
        delbarP(HY_DENS:hyEndPrimVar) = &
             reigP(HY_DENS:hyEndPrimVar,HY_FASTLEFT)*delbarP(HY_FASTLEFT)+&
#ifdef FLASH_USM_MHD
             reigP(HY_DENS:hyEndPrimVar,HY_ALFNLEFT)*delbarP(HY_ALFNLEFT)+&
#endif
             reigP(HY_DENS:hyEndPrimVar,HY_SLOWLEFT)*delbarP(HY_SLOWLEFT)+&
             reigP(HY_DENS:hyEndPrimVar,HY_ENTROPY )*delbarP(HY_ENTROPY )+&
             reigP(HY_DENS:hyEndPrimVar,HY_SLOWRGHT)*delbarP(HY_SLOWRGHT)+&
#ifdef FLASH_USM_MHD
             reigP(HY_DENS:hyEndPrimVar,HY_ALFNRGHT)*delbarP(HY_ALFNRGHT)+&
#endif
             reigP(HY_DENS:hyEndPrimVar,HY_FASTRGHT)*delbarP(HY_FASTRGHT)

        ! delbar at i-1 -------------------------------------------------------------------------
        delbarN(HY_DENS:hyEndPrimVar) = &
             reigN(HY_DENS:hyEndPrimVar,HY_FASTLEFT)*delbarN(HY_FASTLEFT)+&
#ifdef FLASH_USM_MHD
             reigN(HY_DENS:hyEndPrimVar,HY_ALFNLEFT)*delbarN(HY_ALFNLEFT)+&
#endif
             reigN(HY_DENS:hyEndPrimVar,HY_SLOWLEFT)*delbarN(HY_SLOWLEFT)+&
             reigN(HY_DENS:hyEndPrimVar,HY_ENTROPY )*delbarN(HY_ENTROPY )+&
             reigN(HY_DENS:hyEndPrimVar,HY_SLOWRGHT)*delbarN(HY_SLOWRGHT)+&
#ifdef FLASH_USM_MHD
             reigN(HY_DENS:hyEndPrimVar,HY_ALFNRGHT)*delbarN(HY_ALFNRGHT)+&
#endif
             reigN(HY_DENS:hyEndPrimVar,HY_FASTRGHT)*delbarN(HY_FASTRGHT)
     endif

     !! -------------------------------------------------------------------------------------!
     !! [3] Begin polynomial interpolation for PPM ------------------------------------------!
     !! -------------------------------------------------------------------------------------!
     !! initialize arrays with zero
     vecL=0.;    vecR=0.
     !! (a) Parabolic interpolation at the left and right cell interfaces
     !! Colella-Woodward Eqn 1.9, Sekora-Colella Eqn 7, Stone et al Eqn 46
     vecL(HY_DENS:hyEndVar) = 0.5*(Vc(HY_DENS:hyEndVar)+Vm(HY_DENS:hyEndVar)) &
          - (delbar0(HY_DENS:hyEndVar)-delbarN(HY_DENS:hyEndVar))/6.
     vecR(HY_DENS:hyEndVar) = 0.5*(Vc(HY_DENS:hyEndVar)+Vp(HY_DENS:hyEndVar)) &
          - (delbarP(HY_DENS:hyEndVar)-delbar0(HY_DENS:hyEndVar))/6.

     !! (b) PPM interpolate species and mass scalars
#if (NSPECIES+NMASS_SCALARS) > 0
     if (hy_fullSpecMsFluxHandling) then
        vm2Ptr => Smm(:)
        vm1Ptr => Sm (:)
        vc0Ptr => Sc (:)
        vp1Ptr => Sp (:)
        vp2Ptr => Spp(:)

        do isph= 1, HY_NSPEC
           delbarSpP = minmod(vp1Ptr(isph)-vc0Ptr(isph),vp2Ptr(isph)-vp1Ptr(isph))
           delbarSp0 = minmod(vc0Ptr(isph)-vm1Ptr(isph),vp1Ptr(isph)-vc0Ptr(isph))
           delbarSpN = minmod(vm1Ptr(isph)-vm2Ptr(isph),vc0Ptr(isph)-vm1Ptr(isph))

           Sl(isph) = 0.5*(Sc(isph)+vm1Ptr(isph))-(delbarSp0-delbarSpN)/6.
           Sr(isph) = 0.5*(Sc(isph)+vp1Ptr(isph))-(delbarSpP-delbarSp0)/6.

           Sl(isph) = max(min(Sc(isph),vm1Ptr(isph)),min(max(Sc(isph),vm1Ptr(isph)),Sl(isph)))
           Sr(isph) = max(min(Sc(isph),vp1Ptr(isph)),min(max(Sc(isph),vp1Ptr(isph)),Sr(isph)))

           if ((Sr(isph)-Sc(isph))*(Sc(isph)-Sl(isph))<=0.) then
              Sr(isph) = Sc(isph)
              Sl(isph) = Sc(isph)
           endif
           if (6.*(Sr(isph)-Sl(isph))*(Sc(isph)-0.5*(Sr(isph)+Sl(isph))) &
                > (Sr(isph)-Sl(isph))**2) then
              Sl(isph) = 3.*Sc(isph) - 2.*Sr(isph)
           endif
           if (6.*(Sr(isph)-Sl(isph))*(Sc(isph)-0.5*(Sr(isph)+Sl(isph))) &
                <-(Sr(isph)-Sl(isph))**2) then
              Sr(isph) = 3.*Sc(isph) - 2.*Sl(isph)
           endif
           delbarSp0 = Sr(isph)-Sl(isph)
           Sp6 = 6.*(Sc(isph)-0.5*(Sr(isph)+Sl(isph)))

           Sl(isph) = Sl(isph)-min(lambda0(HY_ENTROPY),0.)*hdtn*&
                (delbarSp0+(1.+min(lambda0(HY_ENTROPY),0.)*4./3.*hdtn)*Sp6)

           Sr(isph) = Sr(isph)-max(lambda0(HY_ENTROPY),0.)*hdtn*&
                (delbarSp0-(1.-max(lambda0(HY_ENTROPY),0.)*4./3.*hdtn)*Sp6)
        enddo

        if (hy_flattening) then
           Sl(:) = Flattening*Sc(:) + (1.0-Flattening)*Sl(:)
           Sr(:) = Flattening*Sc(:) + (1.0-Flattening)*Sr(:)
        endif

     endif ! (hy_fullSpecMsFluxHandling)
#endif /*  (NSPECIES+NMASS_SCALARS) > 0 */
     !! End of polynomial interpolation for PPM


     !! -------------------------------------------------------------------------------------!
     !! [4] Contact steepening for PPM  -----------------------------------------------------!
     !! -------------------------------------------------------------------------------------!
     if (hy_ContactSteepening) then
        temp1 = Vp(HY_DENS) - Vm(HY_DENS)
        if (abs(temp1) > hy_tiny) then
           ! Eqn 1.17 : Second derivatives
           del2rhoR = (Vpp(HY_DENS)-2.*Vp(HY_DENS)+ Vc(HY_DENS))/(6.*delta*delta)
           del2rhoL = ( Vc(HY_DENS)-2.*Vm(HY_DENS)+Vmm(HY_DENS))/(6.*delta*delta)

           ! Third derivative
           eta_steep = (del2rhoL-del2rhoR)*delta**2/temp1
           if (del2rhoR*del2rhoL >= 0.) then
              eta_steep = 0.
           endif
           if (epsln*min(Vp(HY_DENS),Vm(HY_DENS))-abs(Vp(HY_DENS) - Vm(HY_DENS)) >= 0.) then
              eta_steep = 0.
           endif

           ! Eqn 1.16
           eta_steep = max(0., min(1., eta1*(eta_steep - eta2)))

           ! Eqn 3.2
           temp2 = abs(Vp(HY_PRES)-Vm(HY_PRES))/min(Vp(HY_PRES),Vm(HY_PRES))
           temp3 = abs(Vp(HY_DENS)-Vm(HY_DENS))/min(Vp(HY_DENS),Vm(HY_DENS))

           if (Vc(HY_GAME)*K0*temp3-temp2 < 0.0) then
              eta_steep = 0.
           endif

           ! Eqn 1.15
           vecL(HY_DENS) = vecL(HY_DENS)*(1.-eta_steep) + (Vm(HY_DENS)+0.5*delbarN(HY_DENS))*eta_steep
           vecR(HY_DENS) = vecR(HY_DENS)*(1.-eta_steep) + (Vp(HY_DENS)-0.5*delbarP(HY_DENS))*eta_steep
        endif
     endif
     !! End of Contact steepening for PPM

     !! -------------------------------------------------------------------------------------!
     !! [5] Flattening for PPM  -------------------------------------------------------------!
     !! -------------------------------------------------------------------------------------!
     if (hy_flattening) then
        vecL(:) = Flattening*Vc(:) + (1.0-Flattening)*vecL(:)
        vecR(:) = Flattening*Vc(:) + (1.0-Flattening)*vecR(:)
     endif

     !! -------------------------------------------------------------------------------------!
     !! [6] Monotonicity check for PPM  -----------------------------------------------------!
     !! -------------------------------------------------------------------------------------!
     ! Ensure that the interpolated values lie between the cell-centered values
     ! Limit according to Colella-Woodward Eqn 1.10
     do n=HY_DENS,hyEndVar
        if ( (vecR(n) - Vc(n))*(Vc(n)-vecL(n)) <= 0.) then
           vecL(n) = Vc(n)
           vecR(n) = Vc(n)
        endif
        if ( 6.*(vecR(n)-vecL(n))*(Vc(n)-0.5*(vecL(n)+vecR(n))) > (vecR(n) - vecL(n))**2  ) then
           vecL(n) = 3.*Vc(n) - 2.*vecR(n)
        endif
        if ( 6.*(vecR(n)-vecL(n))*(Vc(n)-0.5*(vecL(n)+vecR(n))) < -(vecR(n) - vecL(n))**2  ) then
           vecR(n) = 3.*Vc(n) - 2.*vecL(n)
        endif
     enddo
     !! End of Contact steepeing, Flattening, and Monotonicity constraint for PPM


     !! -------------------------------------------------------------------------------------!
     !! [7] Take initial guesses for the left and right states-------------------------------!
     !! -------------------------------------------------------------------------------------!
     !! PPM coefficients for parabolic interpolations
     delW(HY_DENS:hyEndVar) = vecR(HY_DENS:hyEndVar)-vecL(HY_DENS:hyEndVar)
     W6(HY_DENS:hyEndVar)   = 6.*(Vc(HY_DENS:hyEndVar)&
                             -0.5*(vecR(HY_DENS:hyEndVar)+vecL(HY_DENS:hyEndVar)))

     !! [7-a] Right states
     !! Primary variables first
     lambdaMax =max(lambda0(HY_FASTRGHT),0.)
     Wp(HY_DENS:hyEndPrimVar) = vecR(HY_DENS:hyEndPrimVar) - lambdaMax*hdtn &
          *(delW(HY_DENS:hyEndPrimVar) - (1.0 - lambdaMax*hdtn*4./3.)*W6(HY_DENS:hyEndPrimVar))

     !! Secondary variables, gamc, game, eint, 3T vars, grav
     lambdaMax =max(lambda0(HY_ENTROPY),0.)
     Wp(hyEndPrimVar+1:hyEndVar) = vecR(hyEndPrimVar+1:hyEndVar) - lambdaMax*hdtn &
          *(delW(hyEndPrimVar+1:hyEndVar) - (1.0 - lambdaMax*hdtn*4./3.)*W6(hyEndPrimVar+1:hyEndVar))

     !! [7-b] Left states
     !! Primary variables first
     lambdaMin = -min(lambda0(HY_FASTLEFT),0.)
     Wm(HY_DENS:hyEndPrimVar) = vecL(HY_DENS:hyEndPrimVar) + lambdaMin*hdtn &
          *(delW(HY_DENS:hyEndPrimVar) + (1.0 - lambdaMin*hdtn*4./3.)*W6(HY_DENS:hyEndPrimVar))

     !! Secondary variables, gamc, game, eint, 3T vars, grav
     lambdaMin = -min(lambda0(HY_ENTROPY),0.)
     Wm(hyEndPrimVar+1:hyEndVar) = vecL(hyEndPrimVar+1:hyEndVar) + lambdaMin*hdtn &
          *(delW(hyEndPrimVar+1:hyEndVar) + (1.0 - lambdaMin*hdtn*4./3.)*W6(hyEndPrimVar+1:hyEndVar))

     !! [7-c] Apply constraints
     !! Force constant state if simple gamma laws
#ifdef FLASH_EOS_GAMMA
     Wp(HY_GAMC:HY_GAME)=Vc(HY_GAMC:HY_GAME)
     Wm(HY_GAMC:HY_GAME)=Vc(HY_GAMC:HY_GAME)
#endif

     !! [7-d] Gravity component shoud only be spatially reconstructed, 
     !!       without being characteristically traced in time along with fluid's wave
#ifdef GRAVITY
     Wp(HY_GRAV) = vecR(HY_GRAV)
     Wm(HY_GRAV) = vecL(HY_GRAV)
#endif
     !! End of initial guesses
     !! End of high-order polynomial interpolations for PPM interface values


     !! -------------------------------------------------------------------------------------!
     !! [8] Advance the above interpolated interface values by 1/2 time step using ----------!
     !!     characteristic tracing method     -----------------------------------------------!
     !! -------------------------------------------------------------------------------------!
     !! initialize arrays with zero
     sigL=0.;    sigR=0.
     do n=1,HY_WAVENUM
        
        constA = dot_product(leig0(HY_DENS:hyEndPrimVar,n), delW(HY_DENS:hyEndPrimVar))
        constB = dot_product(leig0(HY_DENS:hyEndPrimVar,n),  -W6(HY_DENS:hyEndPrimVar))

        if (hy_RiemannSolver == ROE) then
           if (lambda0(n) < 0.) then
              ! PPM step 10
              !! Left states:
              vecL(HY_DENS:hyEndPrimVar) =  &
                .5*(-1.-   dtn*lambda0(n)                          )*reig0(HY_DENS:hyEndPrimVar,n)*constA &
              +.25*( 1.+2.*dtn*lambda0(n)+4./3.*(dtn*lambda0(n))**2)*reig0(HY_DENS:hyEndPrimVar,n)*constB

              sigL(HY_DENS:hyEndPrimVar) = sigL(HY_DENS:hyEndPrimVar) + vecL(HY_DENS:hyEndPrimVar)

           elseif (lambda0(n) > 0.) then
              !! Right states:
              vecR(HY_DENS:hyEndPrimVar) =  &
                .5*( 1.-   dtn*lambda0(n)                          )*reig0(HY_DENS:hyEndPrimVar,n)*constA &
              +.25*( 1.-2.*dtn*lambda0(n)+4./3.*(dtn*lambda0(n))**2)*reig0(HY_DENS:hyEndPrimVar,n)*constB

              sigR(HY_DENS:hyEndPrimVar) = sigR(HY_DENS:hyEndPrimVar) + vecR(HY_DENS:hyEndPrimVar)

           endif
        else
           !! Left and right states for HLL* type solvers:
           !! For more detail, see "Athena: A new code for astrophysical MHD"
           !! by Stone, Gardiner, Teuben, Hawley, Simon, arXiv:0804.0402v1 [astro-ph] 2 Apr 2008
           !! Apply monotone slope limiting for normal flux
           !! PPM step 10

           !! Left states:
           vecL(HY_DENS:hyEndPrimVar) =  &
                .5*(-1.-   dtn*lambda0(n)                          )*reig0(HY_DENS:hyEndPrimVar,n)*constA &
              +.25*( 1.+2.*dtn*lambda0(n)+4./3.*(dtn*lambda0(n))**2)*reig0(HY_DENS:hyEndPrimVar,n)*constB

           sigL(HY_DENS:hyEndPrimVar) = sigL(HY_DENS:hyEndPrimVar) + vecL(HY_DENS:hyEndPrimVar)

           !! Right states:
           vecR(HY_DENS:hyEndPrimVar) =  &
                .5*( 1.-   dtn*lambda0(n)                          )*reig0(HY_DENS:hyEndPrimVar,n)*constA &
              +.25*( 1.-2.*dtn*lambda0(n)+4./3.*(dtn*lambda0(n))**2)*reig0(HY_DENS:hyEndPrimVar,n)*constB

           sigR(HY_DENS:hyEndPrimVar) = sigR(HY_DENS:hyEndPrimVar) + vecR(HY_DENS:hyEndPrimVar)

        endif
     enddo ! do n=1,HY_WAVENUM

 
     !! -------------------------------------------------------------------------------------!
     !! [9] Consider lower order schemes for 3T variables if requested ----------------------!
     !! -------------------------------------------------------------------------------------!
#ifdef FLASH_UHD_3T
     IF (hy_3Torder .ne. 3) THEN
        if (hy_3Torder == 1) then
           Wp(HY_EINT:HY_ERAD) = Vc(HY_EINT:HY_ERAD)
           Wm(HY_EINT:HY_ERAD) = Vc(HY_EINT:HY_ERAD)

        elseif (hy_3Torder == 2) then


           Wp(HY_EINT:HY_ERAD)=Vc(HY_EINT:HY_ERAD)&
                +0.5*delbar0(HY_EINT:HY_ERAD)*(1.-Flattening)

           Wm(HY_EINT:HY_ERAD)=Vc(HY_EINT:HY_ERAD)&
                -0.5*delbar0(HY_EINT:HY_ERAD)*(1.-Flattening)

           do nVar=HY_EINT,HY_ERAD
              Wm(nVar) = max(min(Vc(nVar),Wm(nVar)),min(max(Vc(nVar),Wm(nVar)),Wm(nVar)))
              Wp(nVar) = max(min(Vc(nVar),Wp(nVar)),min(max(Vc(nVar),Wp(nVar)),Wp(nVar)))
           enddo

           Wp(HY_EINT:HY_ERAD)=Wp(HY_EINT:HY_ERAD)&
                -max(lambda0(HY_ENTROPY),0.)*hdtn*delbar0(HY_EINT:HY_ERAD)*(1.-Flattening)

           Wm(HY_EINT:HY_ERAD)=Wm(HY_EINT:HY_ERAD)&
                -min(lambda0(HY_ENTROPY),0.)*hdtn*delbar0(HY_EINT:HY_ERAD)*(1.-Flattening)
        endif
     ENDIF
#endif


     !! -------------------------------------------------------------------------------------!
     !! [10] Finalize the Riemann states in 1D normal direction -----------------------------!
     !! -------------------------------------------------------------------------------------!
     !! Riemann states in normal direction
     !! Note that gamc, game, grav, & 3T variables are treated separately in the above
     Wm(HY_DENS:hyEndPrimVar) = Vc(HY_DENS:hyEndPrimVar)+W6(HY_DENS:hyEndPrimVar)/12.+sigL(HY_DENS:hyEndPrimVar)
     Wp(HY_DENS:hyEndPrimVar) = Vc(HY_DENS:hyEndPrimVar)+W6(HY_DENS:hyEndPrimVar)/12.+sigR(HY_DENS:hyEndPrimVar)

#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
     Wm(HY_DENS:hyEndPrimVar) = Wm(HY_DENS:hyEndPrimVar)-hdtn*aBnormal(HY_DENS:hyEndPrimVar)*dnBnormal
     Wp(HY_DENS:hyEndPrimVar) = Wp(HY_DENS:hyEndPrimVar)-hdtn*aBnormal(HY_DENS:hyEndPrimVar)*dnBnormal
#ifdef FLASH_UGLM_MHD
     Wm(HY_DENS:hyEndPrimVar) = Wm(HY_DENS:hyEndPrimVar)-hdtn*aGLMnormal(HY_DENS:hyEndPrimVar)*dnGLMnormal
     Wp(HY_DENS:hyEndPrimVar) = Wp(HY_DENS:hyEndPrimVar)-hdtn*aGLMnormal(HY_DENS:hyEndPrimVar)*dnGLMnormal
#endif
#endif

!Wm(HY_MAGX+dir-1)=0.
!Wp(HY_MAGX+dir-1)=0.
     !! DEV-DL: SLOW SHOCK EXPERIMENTS IN THE BELOW ============== 
#ifdef DONGWOOK_SLOWSHOCK
print*,'i am not supposed to be here'
        if ( dir==DIR_X  .and. (Wp(HY_VELX)-Wm(HY_VELX) < 0.) ) then !.and. (Vp(HY_VELX)-Vm(HY_VELX) < 0.)) then
!print*,'one'
           if ( &!(Wp(HY_VELX)-Wm(HY_VELX) < 0.) .or. & !.or. &
                !((Wp(HY_PRES) > Wm(HY_PRES)) .and. (Wp(HY_VELX) > 0. .or. Wm(HY_VELX) > 0.)) .or. &
                !((Wp(HY_PRES) < Wm(HY_PRES)) .and. (Wp(HY_VELX) < 0. .or. Wm(HY_VELX) < 0.)) &
                (( min(abs(lambda(HY_FASTRGHT,lm1)),abs(lambda(HY_FASTRGHT,lc0)))/&
                   max(abs(lambda(HY_FASTRGHT,lm1)),abs(lambda(HY_FASTRGHT,lc0))) >0.9  ) .and. &
                (lambda(HY_FASTRGHT,lm1)>0. .and. lambda (HY_FASTLEFT,lc0)<0.)) .or. &
                (( min(abs(lambda(HY_FASTRGHT,lp1)),abs(lambda(HY_FASTRGHT,lc0)))/&
                   max(abs(lambda(HY_FASTRGHT,lp1)),abs(lambda(HY_FASTRGHT,lc0))) >0.9  ) .and. &
                (lambda (HY_FASTRGHT,lc0)>0. .and. lambda(HY_FASTLEFT,lp1)<0.)) &
                ) then
              Wm(HY_VELX) = Vc(HY_VELX)
              Wp(HY_VELX) = Vc(HY_VELX)
!print*,'hihihi'
!!$              Wm(HY_MAGY) = Vc(HY_MAGY)
!!$              Wp(HY_MAGY) = Vc(HY_MAGY)

           endif
        endif
        if ( dir==DIR_Y .and. (Wp(HY_VELY)-Wm(HY_VELY) < 0.) ) then !.and. (Vp(HY_VELY)-Vm(HY_VELY) < 0.)) then
           if ( &!(Wp(HY_VELY)-Wm(HY_VELY) < 0.) .or. & !.or. &
                !((Wp(HY_PRES) > Wm(HY_PRES)) .and. (Wp(HY_VELY) > 0. .or. Wm(HY_VELY) > 0.)) .or. &
                !((Wp(HY_PRES) < Wm(HY_PRES)) .and. (Wp(HY_VELY) < 0. .or. Wm(HY_VELY) < 0.)) &
                (( min(abs(lambda(HY_FASTRGHT,lm1)),abs(lambda(HY_FASTRGHT,lc0)))/&
                   max(abs(lambda(HY_FASTRGHT,lm1)),abs(lambda(HY_FASTRGHT,lc0))) >0.9  ) .and. &
                (lambda(HY_FASTRGHT,lm1)>0. .and. lambda (HY_FASTLEFT,lc0)<0.)) .or. &
                (( min(abs(lambda(HY_FASTRGHT,lp1)),abs(lambda(HY_FASTRGHT,lc0)))/&
                   max(abs(lambda(HY_FASTRGHT,lp1)),abs(lambda(HY_FASTRGHT,lc0))) >0.9  ) .and. &
                (lambda (HY_FASTRGHT,lc0)>0. .and. lambda(HY_FASTLEFT,lp1)<0.)) &
                ) then
              Wm(HY_VELY) = Vc(HY_VELY)
              Wp(HY_VELY) = Vc(HY_VELY)
           endif
        endif
        if ( dir==DIR_Z ) then
           if ((Wp(HY_VELZ)-Wm(HY_VELZ) < 0.) .or. & !.or. &
                ((Wp(HY_PRES) > Wm(HY_PRES)) .and. (Wp(HY_VELZ) > 0. .or. Wm(HY_VELZ) > 0.)) .or. &
                ((Wp(HY_PRES) < Wm(HY_PRES)) .and. (Wp(HY_VELZ) < 0. .or. Wm(HY_VELZ) < 0.)) &
                ) then
              Wm(HY_VELZ) = Vc(HY_VELZ)
              Wp(HY_VELZ) = Vc(HY_VELZ)
           endif
        endif
        !endif
#endif

#ifdef DWLEE
!!$if ( (Wp(HY_VELX)-Wm(HY_VELX) > 0.) .and. &
!!$     (Vp(HY_VELX)-Vc(HY_VELX) > 0.) .and. &
!!$     (Vc(HY_VELX)-Vm(HY_VELX) > 0.)) then
!!$if ( (Wp(HY_VELX)-Wm(HY_VELX) < 0.) .or. &
!!$     (Vp(HY_VELX)-Vc(HY_VELX) < 0.) .or. &
!!$     (Vc(HY_VELX)-Vm(HY_VELX) < 0.)) then
!print*,"lambda=",lambda0(HY_FASTLEFT),lambda0(HY_FASTRGHT)
     if ( dir==DIR_X ) then !.and. (Wp(HY_VELX)-Wm(HY_VELX) < 0.)) then
        if ( (Wp(HY_VELX)-Wm(HY_VELX) < 0.) .or. & !.or. &
           !((Wp(HY_PRES) > Wm(HY_PRES)) .and. (Wp(HY_VELX) > 0. .or. Wm(HY_VELX) > 0.)) .or. &
           !((Wp(HY_PRES) < Wm(HY_PRES)) .and. (Wp(HY_VELX) < 0. .or. Wm(HY_VELX) < 0.)) &
            (lambda(HY_FASTRGHT,lm1)>0. .and. lambda(HY_FASTLEFT,lc0)<0.) .or. &
            (lambda(HY_FASTRGHT,lc0)>0. .and. lambda(HY_FASTLEFT,lp1)<0.) &
           ) then

           Wm(HY_VELX) = Vc(HY_VELX)
           Wp(HY_VELX) = Vc(HY_VELX)
        endif
     endif
     if ( dir==DIR_Y ) then !.and. (Wp(HY_VELY)-Wm(HY_VELY) < 0.) ) then
        if ( (Wp(HY_VELY)-Wm(HY_VELY) < 0.) .or. & !.or. &
           !((Wp(HY_PRES) > Wm(HY_PRES)) .and. (Wp(HY_VELY) > 0. .or. Wm(HY_VELY) > 0.)) .or. &
           !((Wp(HY_PRES) < Wm(HY_PRES)) .and. (Wp(HY_VELY) < 0. .or. Wm(HY_VELY) < 0.)) &
            (lambda(HY_FASTRGHT,lm1)>0. .and. lambda(HY_FASTLEFT,lc0)<0.) .or. &
            (lambda(HY_FASTRGHT,lc0)>0. .and. lambda(HY_FASTLEFT,lp1)<0.) &
           ) then
           Wm(HY_VELY) = Vc(HY_VELY)
           Wp(HY_VELY) = Vc(HY_VELY)
        endif
     endif
     if ( dir==DIR_Z ) then
        if ((Wp(HY_VELZ)-Wm(HY_VELZ) < 0.) .or. & !.or. &
           ((Wp(HY_PRES) > Wm(HY_PRES)) .and. (Wp(HY_VELZ) > 0. .or. Wm(HY_VELZ) > 0.)) .or. &
           ((Wp(HY_PRES) < Wm(HY_PRES)) .and. (Wp(HY_VELZ) < 0. .or. Wm(HY_VELZ) < 0.)) &
           ) then
           Wm(HY_VELZ) = Vc(HY_VELZ)
           Wp(HY_VELZ) = Vc(HY_VELZ)
        endif
     endif

!$if (shockDetect> 0.) then
!!$Wm(HY_VELX:HY_VELZ) = Vc(HY_VELX:HY_VELZ)
!!$Wp(HY_VELX:HY_VELZ) = Vc(HY_VELX:HY_VELZ)
!!$endif

!!$Wm(HY_VELX:HY_PRES) = Vc(HY_VELX:HY_PRES)
!!$Wp(HY_VELX:HY_PRES) = Vc(HY_VELX:HY_PRES)
#endif

     !! End of advancing the interpolated interface values by 1/2 time step for PPM
  ENDIF ! end of IF (.not. TransUpdateOnly) THEN

End Subroutine Hy_uhd_DataReconstructNormalDir_PPM
