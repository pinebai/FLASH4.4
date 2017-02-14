!!****if* source/physics/Hydro/HydroMain/split/PPM/Hydro_data
!!
!! NAME
!!   
!!  Hydro_data 
!!
!!
!! SYNOPSIS
!!
!!  use Hydro_data
!!
!!
!! DESCRIPTION
!!
!!***


module Hydro_data

#include "Flash.h"
#include "constants.h"

  !!*****Runtime parameters*****
  logical, save :: hy_fluxCorrect,       &
                   hy_hybridRiemann,     &
                   hy_useGravity,        &
                   hy_charLimiting

  
  integer, save :: hy_geometry, &
                   hy_irenorm, &
                   hy_eosMode,  &
                   hy_eosModeAfter, &
                   hy_meshMe, hy_meshNumProcs

  integer, save :: hy_ppmEintFluxConstructionMeth = 0
  integer, save :: hy_ppmEnerFluxConstructionMeth = 0
  integer, save :: hy_ppmEintCFluxConstructionMeth = 0
  integer, save :: hy_ppmEnerCFluxConstructionMeth = 0

  integer, save, dimension(MDIM) :: hy_dirGeom

  real, save :: hy_cfl,         &
                hy_dp_sh_md,    &
                hy_cvisc,       &
                hy_epsiln,      &
                hy_eintSwitch, &
                hy_omg1,        &
                hy_omg2,        &
                hy_vgrid

  logical,save :: hy_useHydro, hy_updateHydroFluxes, &
                  hy_useDiffuse

  !!*****End Runtime parameters*****


  !!*****Directly Derived from Runtime parameters*****
  integer, save :: hy_transverseStencilWidth

  !!*****End Directly Derived from Runtime parameters*****

  ! Maybe should become a runtime parameter.
  ! Set to .TRUE. if you want shock information to be saved in
  ! SHOK_VAR and/or SHKS_VAR, presumably for analyzing shock
  ! locations and strengths, even if hybrid_riemann is off;
  ! see hy_ppm_sweep .
  logical, parameter :: hy_alwaysCallDetectShock = .FALSE.

  logical, save :: hy_movingGrid
  
  !!*****Constants database
  real,    save :: hy_pi


  integer, parameter :: hy_numXN     = NSPECIES+NMASS_SCALARS

  !!******** PPM KERNEL DATA STRUCTURES SCALARS************

  !! - in PPMData
  integer, save :: hy_igodu, hy_nriem, hy_iplm
  real, save    :: hy_small, hy_smallu, hy_smallp, &
                   hy_smlrho, hy_smallx, hy_dp_sh, hy_riemanTol
  logical, save :: hy_ppmModifyStates, hy_leveque,  &
                   hy_useCmaAdvection,          &
                   hy_useCmaFlattening,         &
                   hy_useSteepening,            &
                   hy_useCellAreasForFluxes

  integer,parameter :: lowerFace=1,upperFace=2
  integer, dimension(2,MDIM*2) :: neigh
  !!******** PPM KERNEL DATA STRUCTURES ************

#ifdef FIXEDBLOCKSIZE  
  real,save, DIMENSION(MAXCELLS) :: hy_dela  , &
                                    hy_dp    , hy_du    , hy_dut   , hy_dutt  , &
                                    hy_drho  , hy_dgame , hy_dgamc , hy_dgrav , &
                                    hy_p6    , hy_u6    , hy_ut6   , hy_utt6  , &
                                    hy_rho6  , hy_game6 , hy_gamc6 , hy_grav6 , &
                                    hy_pwl, hy_pwr, hy_dpw, hy_pw6l, hy_pw6r, hy_pwcubic,&
                                    hy_gravl , hy_gravr,                  &
                                    hy_clft  , hy_plft  , hy_uttlft,         &
                                    hy_ulft  , hy_vlft  , hy_utlft ,         &
                                    hy_crght , hy_prght , hy_vrght ,         &
                                    hy_urght , hy_utrght, hy_uttrgt,         &
                                    hy_gmelft, hy_gmergt,                 &
                                    hy_gmclft, hy_gmcrgt,                 &
                                    hy_deint, hy_eint6, hy_eiLft, hy_eiRght

!Ensure array is not zero sized.  This avoids the following 
!link error when compiling with OpenMP on BG/P.
!intrfc.F90:378: undefined reference to `$.TPK$.&&N&&hydro_data%hy_dxn'
#if NSPECIES+NMASS_SCALARS > 0
  real, save, DIMENSION(MAXCELLS,hy_numXN) :: hy_dxn, hy_xn6 ,hy_xnlft, hy_xnrght
#else
  real, save, DIMENSION(MAXCELLS,1) :: hy_dxn, hy_xn6 ,hy_xnlft, hy_xnrght
#endif

  real,save, dimension(2,NYB,NZB,MAXBLOCKS) :: hy_xarea,hy_xdtdx,&
       hy_xgrav,hy_xngrav,hy_xfict
  real,save, dimension(NXB,2,NZB,MAXBLOCKS) :: hy_yarea,hy_ydtdy,&
       hy_ygrav,hy_yngrav,hy_yfict
  real,save, dimension(NXB,NYB,2,MAXBLOCKS) :: hy_zarea,hy_zdtdz,&
       hy_zgrav,hy_zngrav,hy_zfict
  real,save, dimension(2,2,NYB,NZB,MAXBLOCKS) :: hy_xareaAtFaces
  real,save, dimension(2,NXB,2,NZB,MAXBLOCKS) :: hy_yareaAtFaces
  real,save, dimension(2,NXB,NYB,2,MAXBLOCKS) :: hy_zareaAtFaces

#else
  integer, save :: numCells, iguard, jguard, kguard
  real,save,allocatable,dimension(:) :: hy_dela, &
                           hy_dp , hy_du , hy_dut, hy_dutt, &
                           hy_drho, hy_dgame, hy_dgamc, hy_dgrav, &
                           hy_p6, hy_u6 , hy_ut6 , hy_utt6 , &
                           hy_rho6 , hy_game6 , hy_gamc6 , hy_grav6 , &
                           hy_pwl, hy_pwr, hy_dpw, hy_pw6l, hy_pw6r, hy_pwcubic,&
                           hy_gravl, hy_gravr,   &
                           hy_clft  , hy_plft  , hy_uttlft,   &
                           hy_ulft  , hy_vlft  , hy_utlft , &
                           hy_crght , hy_prght , hy_vrght ,         &
                           hy_urght , hy_utrght, hy_uttrgt, &
                           hy_gmelft, hy_gmergt, &
                           hy_gmclft, hy_gmcrgt, &
                           hy_deint, hy_eint6, hy_eiLft, hy_eiRght
  
  real, save, allocatable,dimension(:,:) :: hy_dxn,hy_xn6,hy_xnlft, hy_xnrght
  real,save, allocatable, dimension(:,:,:,:) :: hy_xarea,hy_xdtdx,&
       hy_xgrav,hy_xngrav,hy_xfict
  real,save, allocatable, dimension(:,:,:,:) :: hy_yarea,hy_ydtdy,&
       hy_ygrav,hy_yngrav,hy_yfict
  real,save, allocatable, dimension(:,:,:,:) :: hy_zarea,hy_zdtdz,&
       hy_zgrav,hy_zngrav,hy_zfict
  real,save, allocatable, dimension(:,:,:,:,:) :: hy_xareaAtFaces, &
                                                  hy_yareaAtFaces, &
                                                  hy_zareaAtFaces

#endif

  real,save,allocatable,dimension(:) :: hy_pstor

!! Generated for gravity diagnostic output
  real, save, dimension(MDIM) :: hy_gravMass, hy_gravMassXYZ, hy_gravMassZYX,&
                                 hy_gravMassXZY, hy_gravMassYZX,&
                                 hy_gravMassYXZ, hy_gravMassZXY
  logical,save,dimension(NUNK_VARS) :: hy_gcMask
  integer,save :: hy_gcMaskSize=NUNK_VARS

  character(len=6),save :: hy_FluxRepresentation
  integer, parameter :: hy_numPresFluxes = 1
  integer, save, target, dimension(hy_numPresFluxes) :: hy_specialFluxVars
  data hy_specialFluxVars  / 0 /

  logical, save :: hy_threadBlockList = .false.
  logical, save :: hy_threadWithinBlock = .false.

  !$omp threadprivate(hy_dela,hy_dp,hy_du,hy_dut,hy_dutt,&
  !$omp hy_drho,hy_dgame,hy_dgamc,hy_dgrav,hy_p6,hy_u6,hy_ut6,hy_utt6,&
  !$omp hy_rho6,hy_game6,hy_gamc6,hy_grav6,hy_pwl,hy_pwr,hy_dpw,hy_pw6l,&
  !$omp hy_pw6r,hy_pwcubic,hy_gravl,hy_gravr,hy_clft,hy_plft,hy_uttlft,&
  !$omp hy_ulft,hy_vlft,hy_utlft,hy_crght,hy_prght,hy_vrght,hy_urght,&
  !$omp hy_utrght,hy_uttrgt,hy_gmelft,hy_gmergt,hy_gmclft,hy_gmcrgt,&
  !$omp hy_deint,hy_eint6,hy_eiLft,hy_eiRght,&
  !$omp hy_dxn,hy_xn6,hy_xnlft,hy_xnrght,hy_pstor)

end module Hydro_data
