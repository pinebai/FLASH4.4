!!****if* source/physics/Hydro/HydroMain/unsplit_rad/Hydro_data
!!
!! NAME
!!
!!  Hydro_data
!!
!!
!! SYNOPSIS
!!
!!  MODULE Hydro_data()
!!
!!
!! ARGUMENTS
!!
!!
!! DESCRIPTION
!!
!!  This stores data that are specific to the unsplit Hydro & MHD units.
!!
!!***
 
 
Module Hydro_data

#include "Flash.h"
#include "constants.h"
#include "UHD.h"

  logical, save :: hy_fluxCorrect,        hy_charLimiting,     &
                   hy_shockDetectOn,      hy_shockLowerCFL,&
                   hy_fallbackLowerCFL,                        &
                   hy_useDiffuse,         hy_useViscosity,     &
                   hy_useConductivity,    hy_ContactSteepening,&
                   hy_flattening,         hy_entropy,          &
                   hy_upwindTVD,          hy_use_avisc,        &
                   hy_useGravHalfUpdate,  hy_useGravity,       &
                   hy_updateHydroFluxes,  hy_useHydro,         &
                   hy_use3dFullCTU,       hy_addThermalFlux,   &
                   hy_useHybridOrder,     hy_useVaryingCFL,    &
                   hy_fullRiemannStateArrays,      hy_restart,          &
                   hy_useAuxEintEqn,      hy_hydroComputeDtFirstCall = .true.,&
                   hy_conserveAngMom,                          &
                   hy_fullSpecMsFluxHandling, hy_EOSforRiemann
  logical, save :: hy_computeRadFlBefore = .FALSE.
  logical, save :: hy_doUnsplitLoop0

  logical, save :: hy_useParticles,       hy_useCosmology

  integer, save :: hy_geometry, hy_limiter, hy_irenorm
  integer, save :: hy_unsplitEosMode = MODE_DENS_EI
  integer, save :: hy_EosModeAfter = MODE_DENS_EI
  integer, save :: hy_eosModeGc
  integer, save :: hy_RiemannSolver,hy_RiemannSolverLoc, hy_entropyFixMethod
  integer, save :: hy_hydroComputeDtOption, hy_wenoMethod
  integer, parameter :: hy_numXN = NSPECIES+NMASS_SCALARS
  integer, parameter :: hy_cflStencil = 1

  character(len=MAX_STRING_LENGTH), save :: hy_limiter_str,      &
                                            hy_RiemannSolver_str,&
                                            hy_entropyFixMethod_str,&
                                            hy_wenoMethod_str

  real, save :: hy_cfl
  real, save :: hy_cfl_original
  real, save :: hy_cflFallbackFactor
  real, save :: hy_tiny=1.e-32
  real, save :: hy_Rconst
  real, save :: hy_eswitch
  real, save :: hy_smalldens
  real, save :: hy_smallpres
  real, save :: hy_smallE
  real, save :: hy_LimitedSlopeBeta
  real, save :: hy_maxDiff
  real, save :: hy_cvisc
  real, save :: hy_hybridOrderKappa
  real, save :: hy_fPresInMomFlux
  real, save :: hy_fEddiCorrEarly ! Do Eddington-factor correction early (i.e., in prim2flx)?
  real, save :: hy_lam3ScaleFactor

  logical, save :: hy_reduceGcellFills

  ! Storage for timestep calculation
  integer, save, DIMENSION(5) :: hy_dtminloc
  real, save :: hy_dtminCfl
  real, save :: hy_dtmin
  logical,save :: hy_dtminValid = .FALSE.
  real, save :: hy_dt
  real, save :: hy_simTime       = -HUGE(hy_simTime)
  integer, save :: hy_simGeneration = -1

  ! System of units used
  character(4), save :: hy_units

  ! Everybody should know these!
  integer, save :: hy_meshNumProcs, hy_meshMe
  integer, save :: hy_globalComm

  ! Constants for non-dimensionalization
  real, save :: hy_xref, hy_tref, hy_dref
  real, save :: hy_vref, hy_pref, hy_eref
  real, save :: hy_qref, hy_bref, hy_mref
  real, save :: hy_gref, hy_nref, hy_kref

  ! Order of Accuracy
  ! hy_order=1 ==> First order Godunov
  ! hy_order=2 ==> Muscl-Hancock
  ! hy_order=3 ==> PPM
  ! hy_order=5 ==> WENO
  !! For WENO method, one can also specify wenoMethod="WENO5" or "WENOZ"
  integer, save :: hy_order, hy_transOrder, hy_3Torder
  integer, save :: hy_gcMaskSize
  logical,dimension(NUNK_VARS+NDIM*NFACE_VARS),save :: hy_gcMask

  integer, dimension(NFLUXES), save :: hy_fluxCorVars

  integer, save :: hy_gpotVar = 0
  integer, dimension(MDIM), save :: hy_extraAccelVars = (/0,0,0/)
  logical :: hy_okToUpdateGpot, hy_gpotAlreadyUpToDate

  integer, save :: hy_MOM_PHI = -1, hy_MOM_PHI_FLUX = -1
  integer, save :: hy_MOM_THT = -1, hy_MOM_THT_FLUX = -1
  integer, save :: hy_MOM_ZI = -1, hy_MOM_ZI_FLUX = -1
  integer, save :: hy_MAG_PHI = -1, hy_MAG_PHI_FLUX = -1
  integer, save :: hy_MAG_ZI = -1, hy_MAG_ZI_FLUX = -1

!#ifdef FLASH_UHD_3T
  real,    save :: hy_avogadro, hy_qele
  real,    save :: hy_speedOfLight
  real,    save :: hy_boltz

!#endif

  logical, save :: hy_threadBlockList = .false.
  logical, save :: hy_threadWithinBlock = .false.


!! Temporary array to hold Riemann states without using global SCRATCH arrays

  integer, save :: hy_3TMode = HY3T_NONE

  !! Need for GP interpolation
  real, save :: hy_radiusGP
  real, save :: hy_sigmaGP
  integer, save :: hy_counterGP
  !! GP arrays (we need to initialize hy_WpGP, hy_WmGP, hy_RinvGP)  
  real, save, allocatable, dimension(:,:) :: hy_RinvGP
  real, save, allocatable, dimension(:,:) :: hy_WpGP
  real, save, allocatable, dimension(:,:) :: hy_WmGP



  !! ************************************************************!
  !!                Extra parameters for MHD                    *!
  !! ************************************************************!
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
  logical, save :: hy_forceHydroLimit
#ifdef FLASH_USM_MHD
  logical, save :: hy_killdivb,         &
                   hy_killdivb8w,       &
                   hy_E_modification,   &
                   hy_E_upwind,         &
                   hy_energyFixSwitch
  integer, save :: hy_prol_method = UNDEFINED_PROL
  character(len=MAX_STRING_LENGTH), save :: hy_prol_method_str
#endif

#ifdef FLASH_UGLM_MHD
  real,    save :: hy_C_hyp, hy_C_par
#endif


  logical, save :: hy_useResistivity,   & 
                   hy_useMagneticResistivity
  real,    save :: hy_maxMagDiff

  ! For Biermann Battery Terms
  real,    save :: hy_bier1TZ
  real,    save :: hy_bier1TA
  real,    save :: hy_biermannCoef
  logical, save :: hy_useBiermann1T
  logical, save :: hy_useBiermann
  logical, save :: hy_biermannSource
  logical, save :: hy_hallVelocity
  logical, save :: hy_conserveAngField
#endif

End Module Hydro_data
