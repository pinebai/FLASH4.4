!!****if* source/Grid/GridMain/Chombo/AMR/Grid_data
!!
!! NAME
!!  Grid_data
!!
!! SYNOPSIS
!!
!!  use Grid_data
!!
!! DESCRIPTION 
!!  
!!  This includes the global integer identifier for a block, the grid geometry information
!!  
!!  
!! 
!! CREATE AD:04/12/04
!!
!! 
!!   Defining data structures for storing paramesh related infomation.  
!!   including function for updating the grid information
!!
!! MODIFIED AD:05/19/04
!!   
!!***

Module Grid_data

  implicit none

#include "constants.h"
#include "Flash.h"

!! Define the block information

  integer, save :: gr_iguard = NGUARD
  integer, save :: gr_jguard = NGUARD 
  integer, save :: gr_kguard = NGUARD
  integer, save, dimension(MDIM) :: gr_guard

  integer,save,dimension(MDIM)::gr_bndOrder
  integer,save,dimension(UNK_VARS_BEGIN:UNK_VARS_END) :: gr_vars
  integer,save,dimension(UNK_VARS_BEGIN:UNK_VARS_END) :: gr_vartypes
  logical, save :: gr_anyVarToConvert
  logical, save :: gr_justExchangedGC, gr_allPeriodic, gr_isolatedBoundaries
  logical, save :: gr_useParticles,gr_refineOnParticleCount,gr_refineOnPdens
  integer, save :: gr_minParticlesPerBlk, gr_maxParticlesPerBlk
  integer, save :: gr_globalComm, gr_globalMe, gr_globalNumProcs
  integer, save :: gr_meshComm, gr_meshMe, gr_meshNumProcs
  integer, save :: gr_meshAcrossComm, gr_meshAcrossMe, gr_meshAcrossNumProcs


  integer, save :: gr_eosMode
  integer, save :: gr_eosModeInit, gr_eosModeNow
  integer, save :: gr_oneRefLev=1 !! To be used with the multigrid
  integer ,save :: gr_nrefs
  logical ,save :: gr_convertToConsvdForMeshCalls
  logical ,save :: gr_convertToConsvdInMeshInterp
  logical ,save :: gr_earlyBlockDistAdjustment
  logical, save :: gr_monotone
  integer, save :: gr_intpol
  real, save :: gr_smallx

  integer, save :: gr_numRefineVars, gr_numRefineVarsMax
  integer,allocatable,dimension(:) ,save :: gr_refineVars
  real ,save :: gr_imin,gr_imax,gr_jmin,gr_jmax,gr_kmin,gr_kmax
  real, save, dimension(LOW:HIGH,MDIM) :: gr_globalDomain
  integer,save,dimension(2,MDIM) :: gr_domainBC, gr_blkBC
  integer ,save :: gr_geometry
  integer, save, dimension(MDIM) :: gr_dirGeom
  logical, save, dimension(MDIM) :: gr_dirIsAngular
  logical, save :: gr_geometryOverride
  character(len=MAX_STRING_LENGTH) :: gr_str_geometry
  integer ,save, dimension(MAXREFVARS) :: gr_refine_var    
  real,dimension(MAXREFVARS), save::gr_refine_cutoff,&
       gr_derefine_cutoff,gr_refine_filter
  real, save :: gr_smalle,gr_smallrho
  integer ,dimension(MAXBLOCKS), save :: gr_blkList
  integer, allocatable:: gr_mortons(:)
  real, save :: gr_minCellSize
  real, save, dimension(MDIM) :: gr_minCellSizes

  !below values needed to make data structures for IO output
  integer,save :: gr_globalOffset !stores beginning blk offset for a proc
  integer, save :: gr_globalNumBlocks !

  logical, save :: gr_enableMaskedGCFill

#ifdef GRID_WITH_MONOTONIC
  integer, save :: gr_intpolStencilWidth
#else
  ! The following is for Paramesh3f with native interpolation
  integer, parameter :: gr_intpolStencilWidth = 1
#endif

#ifdef FLASH_PARTICLES


  !gr_blkParticleInfo is a multidimensional array holding various info
  !needed for updating the particle refinement
  !gr_blkParticleInfo(1,:) holds cornerID for 1st dim
  !gr_blkParticleInfo(2,:) holds cornerID for 2nd dim
  !gr_blkParticleInfo(3,:) holds cornerID for 3rd dim
  !gr_blkParticleInfo(4,:) keeps track of if an old block has been processed or not


  integer, save :: gr_blkParticleInfo(4,MAXBLOCKS)
  integer, save :: gr_maxParticlesPerProc
#endif

  integer, save :: gr_numDataStruct
  integer, save, dimension(NDATATYPES) :: gr_gridDataStruct,gr_gridDataStructSize

  integer,save :: gr_lrefineDel, gr_maxRefine
  logical,save :: gr_enforceMaxRefinement
  logical,save :: gr_lrefineMaxRedDoByLogR, gr_lrefineMaxRedDoByTime
  real,save :: gr_lrefineMaxRedRadiusSq
  real,save :: gr_lrefineCenterI,gr_lrefineCenterJ,gr_lrefineCenterK
  real,save :: gr_lrefineMaxRedTimeScale, gr_lrefineMaxRedTRef, gr_lrefineMaxRedLogBase
  integer,save :: gr_restrictAllMethod

  integer, save :: lrefine_min, lrefine_max, grid_changed
  integer, save :: gr_maxBlockSize
  integer, dimension(MDIM), save :: gr_gIndexSize
  real, save :: gr_BRMeshRefineFillRatio
  integer, save :: gr_BRMeshRefineBufferSize
  integer, save :: gr_BRMeshRefineBlockFactor
  integer, save :: gr_tagRadius
  integer, save :: gr_verbosity
  logical, save :: gr_useQuadCFInterp
  logical, save :: gr_useFluxCorrect
  logical, save :: gr_restrictBeforeGhostExchange
  integer, save :: gr_refRatio
  logical, save :: gr_scaleFineFluxes

  real, save, dimension(LOW:HIGH,MDIM) :: gr_region

  logical, save :: gr_bcEnableApplyMixedGds

end Module Grid_data
