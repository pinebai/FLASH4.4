!!****if* source/Grid/GridMain/gr_GCScratchData
!!
!! NAME
!!  gr_GCScratchData
!!
!! SYNOPSIS
!!
!!  use gr_GCScratchData
!!
!! DESCRIPTION  
!!  
!!   This module defines data used for storing the guardcell values of 
!!   center/face datastructures that need to persist between timesteps
!!   for unsplit MHD solvers. The primary purpose of this module is to save
!!   memory.
!! 
!!   
!!***

#include "Flash.h"

Module gr_GCScratchData

  real, save, allocatable, target, dimension(:) :: gr_GCCtr
  real, save, allocatable, target, dimension(:) :: gr_GCFx
  real, save, allocatable, target, dimension(:) :: gr_GCFy
  real, save, allocatable, target, dimension(:) :: gr_GCFz

  integer, save, allocatable, target, dimension(:) :: gr_GCCtrIndList
  integer, save, allocatable, target, dimension(:) :: gr_GCFxIndList
  integer, save, allocatable, target, dimension(:) :: gr_GCFyIndList
  integer, save, allocatable, target, dimension(:) :: gr_GCFzIndList

  integer, save, allocatable, target, dimension(:) :: gr_GCCtrBlkList
  integer, save, allocatable, target, dimension(:) :: gr_GCFxBlkList
  integer, save, allocatable, target, dimension(:) :: gr_GCFyBlkList
  integer, save, allocatable, target, dimension(:) :: gr_GCFzBlkList

  integer, save, dimension(NDIM) :: gr_GCCtrGCCnt
  integer, save, dimension(NDIM) :: gr_GCFxGCCnt
  integer, save, dimension(NDIM) :: gr_GCFyGCCnt
  integer, save, dimension(NDIM) :: gr_GCFzGCCnt

  integer, save :: gr_GCCtrIndCnt
  integer, save :: gr_GCFxIndCnt
  integer, save :: gr_GCFyIndCnt
  integer, save :: gr_GCFzIndCnt

  integer, save :: gr_GCCtrBlkCnt
  integer, save :: gr_GCFxBlkCnt
  integer, save :: gr_GCFyBlkCnt
  integer, save :: gr_GCFzBlkCnt
   
  integer, save, allocatable, target, dimension(:) :: gr_GCCtrBlkOffsets
  integer, save, allocatable, target, dimension(:) :: gr_GCFxBlkOffsets
  integer, save, allocatable, target, dimension(:) :: gr_GCFyBlkOffsets
  integer, save, allocatable, target, dimension(:) :: gr_GCFzBlkOffsets

end Module gr_GCScratchData
