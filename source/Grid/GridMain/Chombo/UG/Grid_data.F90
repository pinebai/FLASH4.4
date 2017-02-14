!!****if* source/Grid/GridMain/Chombo/UG/Grid_data
!!
!! NAME
!!  Grid_data
!!
!! SYNOPSIS
!!
!!  use Grid_data
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

  integer, save :: gr_oneRefLev = 1
  !stores the global number of cells in simulation (w/out guardcells) in each dimension
  integer, save, dimension(MDIM) :: gr_gIndexSize
 
  integer, save, dimension(MDIM) :: gr_procGrid !hold the number of processors in each dim
  integer, save, dimension(MDIM) :: gr_me, gr_guard
  integer, save, dimension(MDIM) :: gr_lIndexSize !stores number of cells in local block (no gcells)

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
  real, save, dimension(LOW:HIGH,MDIM) :: gr_globalDomain

  integer ,save :: gr_iGridSize,gr_jGridSize,gr_kGridSize
  character(len=MAX_STRING_LENGTH), save :: gr_str_geometry
  integer, save :: gr_geometry
  integer, save, dimension(MDIM) :: gr_dirGeom
  logical, save, dimension(MDIM) :: gr_dirIsAngular
  logical, save :: gr_geometryOverride
  integer, save :: gr_eosMode
  integer, save :: gr_eosModeInit, gr_eosModeNow
  logical, save :: gr_justExchangedGC, gr_useParticles, gr_isolatedBoundaries

  integer, save :: gr_verbosity;

  !mostly not used, but convienent for debugging with Grid_dump tool
  integer,save,dimension(UNK_VARS_BEGIN:UNK_VARS_END) :: gr_vars 
 

  integer, save :: gr_iguard,gr_jguard,gr_kguard

  logical, save :: gr_compute_grid_size !used in reordered UG
  real, save :: gr_minCellSize
  real, save, dimension(MDIM) :: gr_minCellSizes

  !below values needed to make data structures for IO output
  integer,save :: gr_globalOffset !stores beginning blk offset for a proc
  integer, save :: gr_globalNumBlocks !
  integer, save, allocatable :: gr_nToLeft(:) !array holding local blocks on each proc, used to calc blk offset for IO
  real, save :: gr_smallx
  real, save :: gr_smalle, gr_smallrho

  real, save, dimension(LOW:HIGH,MDIM) :: gr_region

  logical, save :: gr_bcEnableApplyMixedGds

end Module Grid_data
