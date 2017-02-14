!!****if* source/physics/Hydro/HydroMain/unsplit_rad/hy_uhd_DataReconstructNormalDir_MH
!!
!! NAME
!!
!!  hy_uhd_DataReconstructNormalDir_MH
!!
!! SYNOPSIS
!!
!!  call hy_uhd_DataReconstructNormalDir_MH(integer(IN) :: dir,
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
!!  This subroutine provides a second-order spatially accurate piecewise linear method for
!!  reconstruction. 
!!
!! REFERENCES
!!
!!  * Toro, Riemann Solvers and Numerical Methods for Fluid Dynamics, Springer, 1997
!!  * Lee, D., Ph.D. Dissertation, Univ. of MD, 2006
!!  * Lee, D. and Deane, A., "An Unsplit Staggered Mesh Scheme for Multidimensional
!!                            Magnetohydrodynamics", 228 (2009), 952-975, JCP
!!  * Stone, Gardiner, Teuben, Hawley, Simon, "Athena: A new code for astrophysical MHD"
!!    arXiv:0804.0402v1 [astro-ph] 2 Apr 2008
!!  * Colella, 87, 171-200 (1990), JCP
!!
!!***

Subroutine hy_uhd_DataReconstructNormalDir_MH&
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
                                   hy_3Torder,        &
                                   hy_useAuxEintEqn,  &
                                   hy_fullSpecMsFluxHandling
  use hy_uhd_interface,     ONLY : hy_uhd_TVDslope,   &
                                   hy_uhd_upwindTransverseFlux,&
                                   hy_uhd_eigenParameters, &
                                   hy_uhd_eigenValue,      &
                                   hy_uhd_eigenVector
  use hy_uhd_slopeLimiters, ONLY : mc
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

  integer :: n, nVar, hyBegVar,hyEndVar,hyEndPrimVar
  real    :: dtn,hdtn,factor
  real, dimension(HY_VARINUM)    :: sigL,sigR,vecL,vecR
  real, dimension(HY_VARINUMMAX) :: delbar
  real, pointer, dimension(:) :: vm1Ptr,vc0Ptr,vp1Ptr
  real, dimension(HY_VARINUMMAX), target :: Vm,Vc,Vp

  ! EIG SYSTEM
  logical :: cons=.false.
  real    :: cf,uN,cs,ca,as,af,hyp
  real, dimension(MDIM) :: beta

#if (NSPECIES+NMASS_SCALARS) > 0
  ! Species and mass scalars
  integer :: isph, ispu !isph for hydro scope index, ispu for unk index
  real, dimension(HY_NSPEC), target :: Sc,Sp,Sm
  real, dimension(HY_NSPEC) :: delbarSp 
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


  !! ----------------------------------------------------------------------
  !! Note: MH only requires to use 3-stencil:
  !! 1D array under consideration
  !! ----------------------------------------------------------------------

  !! [A] Usual hydro variables
  !Always initialize with zero before storing arrays!
  Vm=0.;Vc=0.;Vp=0.
  Vm (HY_DENS)         = Data1D(DENS_VAR,         1)
  Vm (HY_VELX:HY_VELZ) = Data1D(VELX_VAR:VELZ_VAR,1)
  Vm (HY_PRES)         = Data1D(PRES_VAR,         1)
  Vm (HY_GAMC:HY_GAME) = Data1D(GAMC_VAR:GAME_VAR,1)
  Vm (HY_EINT)         = Data1D(EINT_VAR,         1)

  Vc (HY_DENS)         = Data1D(DENS_VAR,         2)
  Vc (HY_VELX:HY_VELZ) = Data1D(VELX_VAR:VELZ_VAR,2)
  Vc (HY_PRES)         = Data1D(PRES_VAR,         2)
  Vc (HY_GAMC:HY_GAME) = Data1D(GAMC_VAR:GAME_VAR,2)
  Vc (HY_EINT)         = Data1D(EINT_VAR,         2)

  Vp (HY_DENS)         = Data1D(DENS_VAR,         3)
  Vp (HY_VELX:HY_VELZ) = Data1D(VELX_VAR:VELZ_VAR,3)
  Vp (HY_PRES)         = Data1D(PRES_VAR,         3)
  Vp (HY_GAMC:HY_GAME) = Data1D(GAMC_VAR:GAME_VAR,3)
  Vp (HY_EINT)         = Data1D(EINT_VAR,         3)



  !! [B] magnetic fields for MHD
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
  Vm (HY_MAGX:HY_MAGZ) = Data1D(MAGX_VAR:MAGZ_VAR,1)
  Vc (HY_MAGX:HY_MAGZ) = Data1D(MAGX_VAR:MAGZ_VAR,2)
  Vp (HY_MAGX:HY_MAGZ) = Data1D(MAGX_VAR:MAGZ_VAR,3)
#ifdef FLASH_UGLM_MHD
  Vm (HY_GLMP) = Data1D(HY_GLMP,1)
  Vc (HY_GLMP) = Data1D(HY_GLMP,2)
  Vp (HY_GLMP) = Data1D(HY_GLMP,3)
#endif
#endif

  !! [C] 3T variables
#ifdef FLASH_UHD_3T
#if defined(FLLM_VAR) && defined(HY_FLXL) && defined(HY_PMAT)
  Vm (HY_PRAD)          = Data1D(ERAD_VAR,  1) * Data1D(DENS_VAR,  1) / 3.0
!!$  Vm (HY_FLXL)          = Data1D(FLLM_VAR,         1)
  Vm (HY_EDDI)          =(Data1D(EDDI_VAR,  1) - Data1D(FLLM_VAR,  1)) * Vm(HY_PRAD)
  Vm (HY_PMAT)          = Data1D(PION_VAR,  1) + Data1D(PELE_VAR,  1)
  Vm (HY_PMAT)          = Vm (HY_PMAT) / Vm (HY_PRES)
  Vm (HY_PRAD)          = Vm (HY_PRAD) / Vm (HY_PRES)

  Vc (HY_PRAD)          = Data1D(ERAD_VAR,  2) * Data1D(DENS_VAR,  2) / 3.0
!!$  Vc (HY_FLXL)          = Data1D(FLLM_VAR,         2)
  Vc (HY_EDDI)          =(Data1D(EDDI_VAR,  2) - Data1D(FLLM_VAR,  2)) * Vc(HY_PRAD)
  Vc (HY_PMAT)          = Data1D(PION_VAR,  2) + Data1D(PELE_VAR,  2)
  Vc (HY_PMAT)          = Vc (HY_PMAT) / Vc (HY_PRES)
  Vc (HY_PRAD)          = Vc (HY_PRAD) / Vc (HY_PRES)

  Vp (HY_PRAD)          = Data1D(ERAD_VAR,  3) * Data1D(DENS_VAR,  3) / 3.0
!!$  Vp (HY_FLXL)          = Data1D(FLLM_VAR,         3)
  Vp (HY_EDDI)          =(Data1D(EDDI_VAR,  3) - Data1D(FLLM_VAR,  3)) * Vp(HY_PRAD)
  Vp (HY_PMAT)          = Data1D(PION_VAR,  3) + Data1D(PELE_VAR,  3)
  Vp (HY_PMAT)          = Vp (HY_PMAT) / Vp (HY_PRES)
  Vp (HY_PRAD)          = Vp (HY_PRAD) / Vp (HY_PRES)

#endif
  Vm (HY_EELE:HY_ERAD) = (/Data1D(EELE_VAR,1),Data1D(EION_VAR,1),Data1D(ERAD_VAR,1)/)
  Vc (HY_EELE:HY_ERAD) = (/Data1D(EELE_VAR,2),Data1D(EION_VAR,2),Data1D(ERAD_VAR,2)/)
  Vp (HY_EELE:HY_ERAD) = (/Data1D(EELE_VAR,3),Data1D(EION_VAR,3),Data1D(ERAD_VAR,3)/)
#endif

  !! [D] Gravity component
#ifdef GRAVITY
  Vm (HY_GRAV) = DataGrav1D(1)
  Vc (HY_GRAV) = DataGrav1D(2)
  Vp (HY_GRAV) = DataGrav1D(3)
#endif

  !! [E] Species
#if (NSPECIES+NMASS_SCALARS) > 0
  if (hy_fullSpecMsFluxHandling) then
     do ispu = SPECIES_BEGIN, MASS_SCALARS_END !SPECIES_END
        isph= ispu-NPROP_VARS
        Sm(isph) = Data1D(ispu,1)
        Sc(isph) = Data1D(ispu,2)
        Sp(isph) = Data1D(ispu,3)
     enddo
  endif
#endif


  !!**********************!!
  !! BEGIN RECONSTRUCTION !!
  !!**********************!!

  !! half dt & delta
  dtn=dt/delta
  hdtn=0.5*dtn

  !! initialize arrays with zero
  sigL=0.; sigR=0.
  vecL=0.; vecR=0.
  Wp=0.;   Wm=0.
  !lambda0=0.; reig0=0.; leig0=0.

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
  !! [1] First we compute transverse flux if we are interested in running                 !
  !!     multi-dimensional cases                                                          !
  !! -------------------------------------------------------------------------------------!
  if (NDIM > 1) then

     vm1Ptr => Vm (:)
     vc0Ptr => Vc (:)
     vp1Ptr => Vp (:)

     call hy_uhd_upwindTransverseFlux&
          (dir,hy_transOrder,vm1Ptr,vc0Ptr,vp1Ptr,lambda0,leig0,reig0,HY_END_VARS,sig)

#if (NSPECIES+NMASS_SCALARS) > 0
     if (hy_fullSpecMsFluxHandling) then
        vm1Ptr => Sm (:)
        vc0Ptr => Sc (:)
        vp1Ptr => Sp (:)

        call hy_uhd_upwindTransverseFlux&
          (dir,hy_transOrder,vm1Ptr,vc0Ptr,vp1Ptr,lambda0,leig0,reig0,HY_NSPEC,SpcSig,speciesScalar=.true.)
     endif
#endif
  endif ! NDIM > 1


  !! -------------------------------------------------------------------------------------!
  !! [2] Apply TVD slope limiter for normal gradients                                     !
  !! -------------------------------------------------------------------------------------!
  IF (.not. TransUpdateOnly) THEN
     !! Original
     call hy_uhd_TVDslope(dir,Vm,Vc,Vp,lambda0,leig0,delbar)

     !! Apply flattening if needed
     if (hy_flattening) then
        delbar = delbar*(1.-FlatCoeff)
     endif


  !! -------------------------------------------------------------------------------------!
  !! [3] Apply TVD slope limiter for gamc, game, grav, and 3T variables                   !
  !! -------------------------------------------------------------------------------------!
#ifdef FLASH_EOS_GAMMA /* In the simple case of constant gammas */
     ! We simply use flat-reconstructions as they are constant gammas
     Wp(HY_GAMC:HY_GAME)=Vc(HY_GAMC:HY_GAME)
     Wm(HY_GAMC:HY_GAME)=Vc(HY_GAMC:HY_GAME)
     hyBegVar = HY_EINT
#else
     hyBegVar = HY_GAMC
#endif

     ! Reconstruct game, gamc, eint, grav, and 3T vars
     ! Note: see how hyBegVar is defined in the above.
     do nVar = hyBegVar,hyEndVar !HY_END_VARS
        ! We spatially reconstruct GRAV, & 3T variables here.
        Wp(nVar)=Vc(nVar)+0.5*delbar(nVar)
        Wm(nVar)=Vc(nVar)-0.5*delbar(nVar)

        ! Ensure that the interpolated values lie between the cell-centered values
        Wm(nVar) = max(min(Vc(nVar),Wm(nVar)),min(max(Vc(nVar),Wm(nVar)),Wm(nVar)))
        Wp(nVar) = max(min(Vc(nVar),Wp(nVar)),min(max(Vc(nVar),Wp(nVar)),Wp(nVar)))

        if (nVar .ne. HY_GRAV) then
           ! Gravity components should not be characteristically traced back in time
           ! along with the fluid velocity!
           Wp(nVar)=Wp(nVar)-max(lambda0(HY_ENTROPY),0.)*hdtn*delbar(nVar)
           Wm(nVar)=Wm(nVar)-min(lambda0(HY_ENTROPY),0.)*hdtn*delbar(nVar)
        endif

     enddo

#ifdef FLASH_UGLM_MHD
     nVar     = HY_GLMP
     Wp(nVar) = Vc(nVar)+0.5*delbar(nVar)
     Wm(nVar) = Vc(nVar)-0.5*delbar(nVar)
     Wm(nVar) = max(min(Vc(nVar),Wm(nVar)),min(max(Vc(nVar),Wm(nVar)),Wm(nVar)))
     Wp(nVar) = max(min(Vc(nVar),Wp(nVar)),min(max(Vc(nVar),Wp(nVar)),Wp(nVar)))
#endif



#if (NSPECIES+NMASS_SCALARS) > 0
     if (hy_fullSpecMsFluxHandling) then
        vm1Ptr => Sm (:)
        vc0Ptr => Sc (:)
        vp1Ptr => Sp (:)

        do ispu = SPECIES_BEGIN, MASS_SCALARS_END !SPECIES_END
           isph= ispu-NPROP_VARS
           delbarSp(isph) = mc(vc0Ptr(isph)-vm1Ptr(isph),vp1Ptr(isph)-vc0Ptr(isph))
        enddo

        Sr = Sc+0.5*delbarSp
        Sr = Sr - max(lambda0(HY_ENTROPY),0.)*hdtn*delbarSp

        Sl = Sc-0.5*delbarSp
        Sl = Sl - min(lambda0(HY_ENTROPY),0.)*hdtn*delbarSp

        Sr = max(min(Sc,Sr),min(max(Sc,Sr),Sr))
        Sl = max(min(Sc,Sl),min(max(Sc,Sl),Sl))

     endif ! (hy_fullSpecMsFluxHandling)
#endif /*  (NSPECIES+NMASS_SCALARS) > 0 */

     !! -------------------------------------------------------------------------------------!
     !! [4] Now we perform characteristic tracing to advance by n+1/2                        !
     !! -------------------------------------------------------------------------------------!
     do n=1,HY_WAVENUM
        if (hy_charLimiting) then
           !! Apply slope limiter on characteristic variables

           if (hy_RiemannSolver == ROE) then
              if (lambda0(n) < 0.) then

                 !! Apply monotone slope limiting for normal flux
                 vecL(HY_DENS:hyEndPrimVar) = .5*(-1.-dtn*lambda0(n))*reig0(HY_DENS:hyEndPrimVar,n)*delbar(n)
                 sigL(HY_DENS:hyEndPrimVar) = sigL(HY_DENS:hyEndPrimVar) + vecL(HY_DENS:hyEndPrimVar)
              else
                 !! Apply monotone slope limiting for normal flux term
                 vecR(HY_DENS:hyEndPrimVar) = .5*(1.-dtn*lambda0(n))*reig0(HY_DENS:hyEndPrimVar,n)*delbar(n)
                 sigR(HY_DENS:hyEndPrimVar) = sigR(HY_DENS:hyEndPrimVar) + vecR(HY_DENS:hyEndPrimVar)
              endif
           else
              !! Left states for HLL* type solvers
              !! For more detail, see "Athena: A new code for astrophysical MHD"
              !! by Stone, Gardiner, Teuben, Hawley, Simon, arXiv:0804.0402v1 [astro-ph] 2 Apr 2008
              !! Apply monotone slope limiting for normal flux
              vecL(HY_DENS:hyEndPrimVar) = .5*(-1.-dtn*lambda0(n))*reig0(HY_DENS:hyEndPrimVar,n)*delbar(n)
              sigL(HY_DENS:hyEndPrimVar) = sigL(HY_DENS:hyEndPrimVar) + vecL(HY_DENS:hyEndPrimVar)

              !! Right states for HLL* type solvers
              !! Apply monotone slope limiting for normal flux term
              vecR(HY_DENS:hyEndPrimVar) = .5*(1.-dtn*lambda0(n))*reig0(HY_DENS:hyEndPrimVar,n)*delbar(n)
              sigR(HY_DENS:hyEndPrimVar) = sigR(HY_DENS:hyEndPrimVar) + vecR(HY_DENS:hyEndPrimVar)


#ifdef FLASH_UGLM_MHD
!!$              sigL(HY_GLMP) = 0.
!!$              sigR(HY_GLMP) = 0.

!!$              vecL(HY_GLMP) = -0.5*reig0(HY_GLMP,n)*delbar(n)
!!$              vecR(HY_GLMP) =  0.5*reig0(HY_GLMP,n)*delbar(n)
!!$              sigL(HY_GLMP) =  sigL(HY_GLMP)+vecL(HY_GLMP)
!!$              sigR(HY_GLMP) =  sigR(HY_GLMP)+vecR(HY_GLMP)

!!$              if (dir == DIR_X) then
!!$                 sigL(HY_MAGX) = 0.
!!$                 sigR(HY_MAGX) = 0.
!!$
!!$                 vecL(HY_MAGX) = -0.5*reig0(HY_MAGX,n)*delbar(HY_MAGX)
!!$                 vecR(HY_MAGX) =  0.5*reig0(HY_MAGX,n)*delbar(HY_MAGX)
!!$                 sigL(HY_MAGX) =  sigL(HY_MAGX)+vecL(HY_MAGX)
!!$                 sigR(HY_MAGX) =  sigR(HY_MAGX)+vecR(HY_MAGX)
!!$              elseif (dir == DIR_Y) then
!!$                 sigL(HY_MAGY) = 0.
!!$                 sigR(HY_MAGY) = 0.
!!$
!!$                 vecL(HY_MAGY) = -0.5*reig0(HY_MAGY,n)*delbar(HY_MAGY)
!!$                 vecR(HY_MAGY) =  0.5*reig0(HY_MAGY,n)*delbar(HY_MAGY)
!!$                 sigL(HY_MAGY) =  sigL(HY_MAGY)+vecL(HY_MAGY)
!!$                 sigR(HY_MAGY) =  sigR(HY_MAGY)+vecR(HY_MAGY)
!!$              endif

#endif

           endif

        else
           !! Apply slope limiter on primitive variables
           if (hy_RiemannSolver == ROE) then
              if (lambda0(n) < 0.) then
                 !! Apply monotone slope limiting for normal flux
                 vecL(HY_DENS:hyEndPrimVar) = .5*(-1.- dtn*lambda0(n))*reig0(HY_DENS:hyEndPrimVar,n)&
                      *dot_product(leig0(HY_DENS:hyEndPrimVar,n), delbar(HY_DENS:hyEndPrimVar))
                 sigL(HY_DENS:hyEndPrimVar) = sigL(HY_DENS:hyEndPrimVar) + vecL(HY_DENS:hyEndPrimVar)

              else
                 !! Apply monotone slope limiting for normal flux
                 vecR(HY_DENS:hyEndPrimVar) = .5*(1.-dtn*lambda0(n))*reig0(HY_DENS:hyEndPrimVar,n)&
                      *dot_product(leig0(HY_DENS:hyEndPrimVar,n), delbar(HY_DENS:hyEndPrimVar))
                 sigR(HY_DENS:hyEndPrimVar) = sigR(HY_DENS:hyEndPrimVar) + vecR(HY_DENS:hyEndPrimVar)
              endif
           else
              !! Left states for HLL* type solvers
              !! Apply monotone slope limiting for normal flux
              vecL(HY_DENS:hyEndPrimVar) = .5*(-1.- dtn*lambda0(n))*reig0(HY_DENS:hyEndPrimVar,n)&
                   *dot_product(leig0(HY_DENS:hyEndPrimVar,n), delbar(HY_DENS:hyEndPrimVar))
              sigL(HY_DENS:hyEndPrimVar) = sigL(HY_DENS:hyEndPrimVar) + vecL(HY_DENS:hyEndPrimVar)

              !! Right states for HLL* type solvers
              !! Apply monotone slope limiting for normal flux
              vecR(HY_DENS:hyEndPrimVar) = .5*(1.-dtn*lambda0(n))*reig0(HY_DENS:hyEndPrimVar,n)&
                   *dot_product(leig0(HY_DENS:hyEndPrimVar,n), delbar(HY_DENS:hyEndPrimVar))
              sigR(HY_DENS:hyEndPrimVar) = sigR(HY_DENS:hyEndPrimVar) + vecR(HY_DENS:hyEndPrimVar)

#ifdef FLASH_UGLM_MHD
!!$              sigL(HY_GLMP) = 0.
!!$              sigR(HY_GLMP) = 0.

!!$              vecL(HY_GLMP) = -0.5*reig0(HY_GLMP,n)*leig0(HY_GLMP,n)*delbar(HY_GLMP)
!!$              vecR(HY_GLMP) =  0.5*reig0(HY_GLMP,n)*leig0(HY_GLMP,n)*delbar(HY_GLMP)
!!$              sigL(HY_GLMP) =  sigL(HY_GLMP)+vecL(HY_GLMP)
!!$              sigR(HY_GLMP) =  sigR(HY_GLMP)+vecR(HY_GLMP)

!!$              if (dir == DIR_X) then
!!$                 sigL(HY_MAGX) = 0.
!!$                 sigR(HY_MAGX) = 0.
!!$
!!$                 vecL(HY_MAGX) = -0.5*reig0(HY_MAGX,n)*leig0(HY_MAGX,n)*delbar(HY_MAGX)
!!$                 vecR(HY_MAGX) =  0.5*reig0(HY_MAGX,n)*leig0(HY_MAGX,n)*delbar(HY_MAGX)
!!$                 sigL(HY_MAGX) =  sigL(HY_MAGX)+vecL(HY_MAGX)
!!$                 sigR(HY_MAGX) =  sigR(HY_MAGX)+vecR(HY_MAGX)
!!$              elseif (dir == DIR_Y) then
!!$                 sigL(HY_MAGY) = 0.
!!$                 sigR(HY_MAGY) = 0.
!!$
!!$                 vecL(HY_MAGY) = -0.5*reig0(HY_MAGY,n)*leig0(HY_MAGY,n)*delbar(HY_MAGY)
!!$                 vecR(HY_MAGY) =  0.5*reig0(HY_MAGY,n)*leig0(HY_MAGY,n)*delbar(HY_MAGY)
!!$                 sigL(HY_MAGY) =  sigL(HY_MAGY)+vecL(HY_MAGY)
!!$                 sigR(HY_MAGY) =  sigR(HY_MAGY)+vecR(HY_MAGY)
!!$              endif
#endif


           endif

        endif ! End of if (hy_charLimiting)
     enddo ! End of do n=1,HY_WAVENUM


     !! -------------------------------------------------------------------------------------!
     !! [5] Apply first-order if requested                                                   !
     !! -------------------------------------------------------------------------------------!
#ifdef FLASH_UHD_3T
     if (hy_3Torder == 1) then
        Wp(HY_EINT:HY_ERAD) = Vc(HY_EINT:HY_ERAD)
        Wm(HY_EINT:HY_ERAD) = Vc(HY_EINT:HY_ERAD)
     endif
#endif


     !! -------------------------------------------------------------------------------------!
     !! [6] Finalize the Riemann states in 1D normal direction                               !
     !! -------------------------------------------------------------------------------------!
     !! Riemann states in normal direction
     !! Note that gamc, game, grav, & 3T variables are treated separately in the above
     Wm(HY_DENS:hyEndPrimVar) = Vc(HY_DENS:hyEndPrimVar)+sigL(HY_DENS:hyEndPrimVar)
     Wp(HY_DENS:hyEndPrimVar) = Vc(HY_DENS:hyEndPrimVar)+sigR(HY_DENS:hyEndPrimVar)

#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
     Wm(HY_DENS:hyEndPrimVar) = Wm(HY_DENS:hyEndPrimVar)-hdtn*aBnormal(HY_DENS:hyEndPrimVar)*dnBnormal
     Wp(HY_DENS:hyEndPrimVar) = Wp(HY_DENS:hyEndPrimVar)-hdtn*aBnormal(HY_DENS:hyEndPrimVar)*dnBnormal
#ifdef FLASH_UGLM_MHD
     Wm(HY_DENS:hyEndPrimVar) = Wm(HY_DENS:hyEndPrimVar)-hdtn*aGLMnormal(HY_DENS:hyEndPrimVar)*dnGLMnormal
     Wp(HY_DENS:hyEndPrimVar) = Wp(HY_DENS:hyEndPrimVar)-hdtn*aGLMnormal(HY_DENS:hyEndPrimVar)*dnGLMnormal
#endif
#endif

!!$
!!$     Wm(HY_GLMP) = Vc(HY_GLMP)
!!$     Wp(HY_GLMP) = Vc(HY_GLMP)


  ENDIF ! end of IF (.not. TransUpdateOnly) THEN




End Subroutine Hy_uhd_DataReconstructNormalDir_MH
