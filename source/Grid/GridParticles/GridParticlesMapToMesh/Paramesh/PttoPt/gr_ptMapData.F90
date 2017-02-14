!!****if* source/Grid/GridParticles/GridParticlesMapToMesh/Paramesh/PttoPt/gr_ptMapData
!!
!! NAME
!!
!!  gr_ptMapData
!!
!! SYNOPSIS
!!
!!  use gr_ptMapData
!!
!! DESCRIPTION
!!
!!  Data module for the variables used to map particles to the grid.
!!
!! ARGUMENTS
!!
!! NOTES 
!!
!!  When compiling on an IBM platform let the compiler perform the appropriate
!!  alignment.  The same is true for any RISC processor.  Do not add the 
!!  "sequence" keyword to the data structure as this will override any 
!!  requested compiler alignment.
!!  We require natural data alignment, i.e 8 byte values on 8 byte boundary.
!!
!!***

#include "constants.h"
#include "Flash.h"

module gr_ptMapData

  implicit none

!gr_ptSmearLen describes how many cells from the central cell containing the 
!particle can possibly receive mass accumulation.
!e.g.  NGP scheme: gr_ptSmearLen = 0, CIC/TSC scheme: gr_ptSmearLen = 1
  integer, save :: gr_ptSmearLen
  integer, save, allocatable, dimension(:) :: gr_ptRecvSpecifier, &
       gr_ptRecvSpecifierTmp, gr_ptRecvTotal, gr_ptRecvTotalTmp
  integer, save :: gr_ptNumMessagesToSend

#ifdef FLASH_GRID_PARAMESH

  integer, parameter :: MAXNEGH = 2**(NDIM-1)  !1(1D), 2(2D), 4(3D)
  integer, parameter :: NUMBGUARDREGIONS = (3**NDIM)-1

  !Hierarchy is: SourceBlock -> GuardCellRegion -> GuardCellRegionNeighbor
  !NOTE: Do not add the "sequence" keyword.
  type GuardCellRegionNeighbor
     integer, dimension(3) :: negh !Encloses blockID, procID, refLevel
     integer, dimension(NDIM) :: cornerID
     integer, dimension(LOW:HIGH, NDIM) :: srcCoords, destCoords
  end type GuardCellRegionNeighbor

  type GuardCellRegion
     type (GuardCellRegionNeighbor), dimension(MAXNEGH) :: neighbor
     integer :: numNegh
  end type GuardCellRegion

  !Each source block contains 2, 8, 26 guard cell regions.
  type SourceBlock
     type (GuardCellRegion), dimension(NUMBGUARDREGIONS) :: haloRegion     
     real, dimension(LOW:HIGH,NDIM) :: bndBlk
     real, dimension(1:NDIM) :: cellSpacing
     integer :: blockID
  end type SourceBlock

  !The domain on each processor will change each time step.
  type (SourceBlock), save, allocatable, dimension(:) :: gr_ptDomain

  !Total amount of memory in 3D: 
  !GuardCellRegionNeighbor = (18 * 4) bytes = 72 bytes
  !GuardCellRegion = ((4 * 72) + 4) bytes = 292 bytes
  !SourceBlock = ((26 * 292) + (6 * 8)) bytes = 7640 bytes
  !gr_ptDomain (say 1000 blocks) = (7640 * 1000) bytes = about 8 Megabytes at best 

#endif

end module gr_ptMapData
