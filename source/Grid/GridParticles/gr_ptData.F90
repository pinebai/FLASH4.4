!!****if* source/Grid/GridParticles/gr_ptData
!!
!! NAME
!!
!!  gr_ptData
!!
!! SYNOPSIS
!!
!!  use gr_ptData
!!
!! DESCRIPTION
!!
!!  Data module for the variables used in the MoveSieve version of
!!  GridParticles implementation for Paramesh
!!
!! ARGUMENTS
!!
!!
!!***
#include "Flash.h"

Module gr_ptData

  implicit none
  integer,save,dimension(MAXBLOCKS) :: gr_ptBlkList
  integer,save :: gr_ptBlkCount,gr_ptRemoveAlgo,gr_ptNumToReduce
  integer,save :: gr_ptMaxPerProc, gr_ptSieveCheckFreq, gr_ptSieveFreq
  integer,save :: gr_ptMaxVirtualCount
  integer,save :: gr_ptLogLevel
  real,save,allocatable,dimension(:,:) :: gr_ptDestBuf,gr_ptSourceBuf
  real,save,allocatable,dimension(:,:,:) :: gr_ptBuf

  logical,save :: gr_ptRemove, gr_ptKeepLostParticles
  
  integer, save :: gr_ptBlk, gr_ptProc, gr_ptTag, &
       gr_ptPosx, gr_ptPosy, gr_ptPosz, &
       gr_ptPos2x, gr_ptPos2y, gr_ptPos2z, &
       gr_ptVelx, gr_ptVely, gr_ptVelz, &
       gr_ptVel2x, gr_ptVel2y, gr_ptVel2z, gr_ptVirtual
  logical, save :: gr_ptPosTmp, gr_ptVel, gr_ptVelTmp

  integer, save :: gr_ptToggle, gr_ptProcDist
  integer, save :: gr_ptSmearLen=1

#ifndef FLASH_GRID_UG
  real, save ::    gr_ptMaxPerProcUpperThresh
  real, save ::    gr_ptMaxPerProcLowerThresh
  real, save ::    gr_ptMaxPerProcBlockFactor
  integer, save :: gr_ptMaxPerProcBlockNoFuzz
  logical, save :: gr_ptRefineOnPtMaxPerProc
#endif  
end Module gr_ptData
