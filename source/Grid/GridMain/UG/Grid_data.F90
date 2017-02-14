!!****if* source/Grid/GridMain/UG/Grid_data
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
!!  Defining data structures for Uniform Grid
!!
!!***

!!REORDER(5):scratch, scratch_ctr, scratch_facevar[xyz]


Module Grid_data
  implicit none

#include "Flash.h"
#include "constants.h"

#ifdef FIXEDBLOCKSIZE


#if NSCRATCH_GRID_VARS > 0
  real, save, target :: scratch(SCRATCH_GRID_VARS_BEGIN:SCRATCH_GRID_VARS_END,&
                                GRID_ILO_GC:GRID_IHI_GC+1, &
                                GRID_JLO_GC:GRID_JHI_GC+1, &
                                GRID_KLO_GC:GRID_KHI_GC+1,1)
       
#else
  real, save,target :: scratch(1,1,1,1,1)
#endif

#if NSCRATCH_CENTER_VARS > 0
  real, save, target :: scratch_ctr(SCRATCH_CENTER_VARS_BEGIN:SCRATCH_CENTER_VARS_END,&
                 GRID_ILO_GC:GRID_IHI_GC,&
                 GRID_JLO_GC:GRID_JHI_GC,&
                 GRID_KLO_GC:GRID_KHI_GC,1)
#else
  real, save,target :: scratch_ctr(1,1,1,1,1)
#endif

#if(NSCRATCH_FACEX_VARS>0)
  real, save, target :: &
       scratch_facevarx( SCRATCH_FACEX_VARS_BEGIN:SCRATCH_FACEX_VARS_END,&
                         GRID_ILO_GC:GRID_IHI_GC+1,&
                         GRID_JLO_GC:GRID_JHI_GC,  &
                         GRID_KLO_GC:GRID_KHI_GC,1)
#else
  real, save,target,dimension(1,1,1,1,1)::scratch_facevarx
#endif

#if(NSCRATCH_FACEY_VARS>0)
  real, save, target :: &
       scratch_facevary( SCRATCH_FACEY_VARS_BEGIN:SCRATCH_FACEY_VARS_END,&
                         GRID_ILO_GC:GRID_IHI_GC,    &
                         GRID_JLO_GC:GRID_JHI_GC+K2D,&
                         GRID_KLO_GC:GRID_KHI_GC,1)
#else
  real, save,target,dimension(1,1,1,1,1)::scratch_facevary
#endif

#if(NSCRATCH_FACEZ_VARS>0)
  real, save, target :: &
       scratch_facevarz( SCRATCH_FACEZ_VARS_BEGIN:SCRATCH_FACEZ_VARS_END,&
                         GRID_ILO_GC:GRID_IHI_GC, &
                         GRID_JLO_GC:GRID_JHI_GC, &
                         GRID_KLO_GC:GRID_KHI_GC+K3D,1) 
#else
  real, save,target,dimension(1,1,1,1,1)::scratch_facevarz
#endif


#else
  !NONFIXEDBLOCKSIZE.
  real, save, target, allocatable :: scratch(:,:,:,:,:)
  real, save, target, allocatable :: scratch_ctr(:,:,:,:,:)
  real, save, target, allocatable :: scratch_facevarx(:,:,:,:,:)
  real, save, target, allocatable :: scratch_facevary(:,:,:,:,:)
  real, save, target, allocatable :: scratch_facevarz(:,:,:,:,:)
#endif
!End of ifdef FIXEDBLOCKSIZE


  integer, parameter :: lnblocks = 1
  integer, parameter :: SAVED_VARS=1
  integer, save :: gr_oneRefLev = 1
  !stores the global number of cells in simulation (w/out guardcells) in each dimension
  integer, save, dimension(MDIM) :: gr_gIndexSize
 
  integer, save, dimension(MDIM) :: gr_guard
  integer, save, dimension(MDIM) :: gr_lIndexSize !stores number of cells in local block (no gcells)

  !stores lower left hand global cell index for each local block in each dim
  integer, save, dimension(MDIM) :: gr_blkCornerID
  integer, save :: gr_blockType

  integer, save :: gr_meshComm, gr_meshMe, gr_meshNumProcs
  integer, save ::  gr_meshAcrossComm, gr_meshAcrossMe, gr_meshAcrossNumProcs
  integer, save :: gr_globalComm, gr_globalMe, gr_globalNumProcs
  integer, save, dimension(MDIM) :: gr_axisComm, gr_axisMe, gr_axisNumProcs
  integer, save, dimension(NDATATYPES,MDIM) :: gr_offset
  integer, save, dimension(NDATATYPES,MDIM) :: gr_exch
  integer, save :: gr_numDataStruct
  integer, save, dimension(NDATATYPES) :: gr_gridDataStruct,gr_gridDataStructSize

  integer, save, dimension(MDIM) :: gr_bndOrder
  logical, save :: gr_allPeriodic
  integer, save, dimension(1)::gr_blkList
  integer, save, dimension(2,MDIM) :: gr_domainBC, gr_blkBC
  real ,save :: gr_imin,gr_imax,gr_jmin,gr_jmax,gr_kmin,gr_kmax

  integer ,save :: gr_iGridSize,gr_jGridSize,gr_kGridSize
  character(len=MAX_STRING_LENGTH), save :: gr_str_geometry
  integer, save :: gr_geometry
  integer, save, dimension(MDIM) :: gr_dirGeom
  logical, save, dimension(MDIM) :: gr_dirIsAngular
  logical, save :: gr_geometryOverride
  integer, save :: gr_eosMode
  integer, save :: gr_eosModeInit, gr_eosModeNow
  logical, save :: gr_justExchangedGC, gr_useParticles, gr_isolatedBoundaries

  logical, save :: gr_useEnergyDeposition

  !mostly not used, but convienent for debugging with Grid_dump tool
  integer,save,dimension(UNK_VARS_BEGIN:UNK_VARS_END) :: gr_vars 
 
  real, save, dimension(MDIM) :: gr_delta
  real, save, dimension(LOW:HIGH,MDIM) :: gr_globalDomain

  integer, save :: gr_iguard,gr_jguard,gr_kguard
  integer, save :: gr_iloGc,gr_ihiGc,gr_jloGc,gr_jhiGc,gr_kloGc,gr_khiGc
  integer, save :: gr_ilo,gr_ihi,gr_jlo,gr_jhi,gr_klo,gr_khi
  real, save,target, allocatable, dimension(:,:,:) ::  gr_iCoords,gr_jCoords,gr_kCoords

  logical, save :: gr_compute_grid_size !used in reordered UG
  real, save :: gr_minCellSize
  real, save, dimension(MDIM) :: gr_minCellSizes

  !below values needed to make data structures for IO output
  integer,save :: gr_globalOffset !stores beginning blk offset for a proc
  integer,save,target,allocatable, dimension(:,:) :: gr_gid  !holds neigh, child, parent info for checkpoint files
  integer, save :: gr_globalNumBlocks !
  integer, save, allocatable :: gr_nToLeft(:) !array holding local blocks on each proc, used to calc blk offset for IO
  real, save :: gr_smallx
  real, save :: gr_smalle, gr_smallrho
  logical, save :: gr_refineOnParticleCount=.false.,gr_refineOnPdens=.false.

  real, save, dimension(LOW:HIGH,MDIM) :: gr_region
  integer, save :: surr_blks(3,3,1+2*K2D,1+2*K3D,1)

  logical, save :: gr_bcEnableApplyMixedGds

end Module Grid_data
