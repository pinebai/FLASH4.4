!!****if* source/Grid/GridSolvers/Pfft/gr_pfftData
!!
!! NAME
!!  gr_pfftData
!!
!! SYNOPSIS
!!
!!  use gr_pfftData
!!
!!  This includes subunit scope data for the pfft solver
!!  
!!
!!***

module gr_pfftData
#include "constants.h"
#include "Flash.h"
  implicit none
  integer, save, dimension(MDIM) :: pfft_globalLen, pfft_inLen, pfft_t1Len, &
       pfft_midLen,pfft_t2Len,pfft_transformType, pfft_comm, &
       pfft_me, pfft_procGrid, pfft_dimOrder
  integer, save, target :: pfft_outLen(MDIM)
  integer, save, target :: pfft_pclbaseDatType(0:MDIM)
  integer, save, target :: pfft_localLimits(LOW:HIGH,MDIM)
  integer, save, dimension(LOW:HIGH,MDIM) :: pfft_pclInGcNumGlobal,&
       pfft_pclMidGcNumGlobal, pfft_pclOutGcNumGlobal
  integer, save, dimension(LOW:HIGH,MDIM) :: pfft_pclInGcNumLocal,&
       pfft_pclMidGcNumLocal, pfft_pclOutGcNumLocal
  real, save, allocatable, dimension(:) :: pfft_work1, pfft_work2
  real, save, pointer, dimension(:) :: pfft_trigIaxis, pfft_trigJaxis, &
       pfft_trigKaxis
  real, save, dimension(MDIM) :: pfft_scale
  integer, save :: pfft_workSize, pfft_ndim, pfft_commWithTopology, & 
       pfft_myPE, pfft_numProcs, pfft_mode
  logical, save :: pfft_needMap,pfft_setupOnce, pfft_usableProc
  real, save, allocatable,dimension(:,:) :: pfft_wave
  integer, save :: gr_pfftDiffOpDiscretize
  logical, save :: pfft_isSimpleNonMappedAMR1D, pfft_inRegion
  real, dimension(LOW:HIGH,MDIM), save :: pfft_regionBndBox
end module gr_pfftData
