!!****if* source/physics/Gravity/GravityMain/PlanePar/Gravity_accelOneRow
!!
!! NAME
!!
!!  Gravity_accelOneRow 
!!
!! SYNOPSIS
!!
!!  call Gravity_accelOneRow(integer(IN)  :: pos(2),
!!                           integer(IN)  :: sweepDir,
!!                           integer(IN)  :: blockID,
!!                           integer(IN)  :: numCells,
!!                           real(INOUT)  :: grav(numCells),
!!                           integer(IN),optional :: potentialIndex,
!!                           integer(IN),optional :: extraAccelVars(MDIM))
!!
!! DESCRIPTION
!!
!!  This routine computes the gravitational acceleration for a row
!!  of cells in a specified direction in a given block. This
!!  version implements a plane-parallel gravitational field.
!!
!! ARGUMENTS
!!
!!  pos      :  Row indices transverse to the sweep direction
!!  sweepDir :    The sweep direction:  allowed values are 
!!              SWEEP_X, SWEEP_Y, and SWEEP_Z. These values are defined
!!              in constants.h.
!!  blockID  :  The local identifier of the block to work on
!!  numCells :  Number of cells to update in grav()
!!  grav()   :   Array to receive result
!!  potentialIndex :  optional, not applicable in planepar gravity
!!  extraAccelVars :  optional, ignored in this implementation
!! 
!!***

subroutine Gravity_accelOneRow (pos, sweepDir, blockID, numCells, grav, &
                                potentialIndex, extraAccelVars)

!========================================================================

  use Gravity_data, ONLY: useGravity
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_getCellCoords

  implicit none

#include "Flash.h"
#include "constants.h"

  integer, intent(IN) :: sweepDir,blockID,numCells
  integer, dimension(2),INTENT(in) ::pos
  real, dimension(numCells),INTENT(inout) :: grav
  integer,intent(IN),optional :: potentialIndex
  integer,intent(IN),OPTIONAL :: extraAccelVars(MDIM)

  real, DIMENSION(MDIM) :: grav_zone
  real, DIMENSION(MDIM) :: gc, gl, gr

#ifdef FIXEDBLOCKSIZE
  real,dimension(GRID_KHI_GC) :: zCenter,zLeft,zRight
  real,dimension(GRID_JHI_GC) :: yCenter,yLeft,yRight
  real,dimension(GRID_IHI_GC) :: xCenter,xLeft,xRight
#else
  real,allocatable,dimension(:) ::xCenter,xLeft,xRight
  real,allocatable,dimension(:) ::yCenter,yLeft,yRight
  real,allocatable,dimension(:) ::zCenter,zLeft,zRight
#endif
  
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
  integer :: i,j,k
  integer :: sizeX,sizeY,sizeZ
  logical :: gcell=.true.

!
!==============================================================================
  if (.NOT.useGravity) return

#ifndef FIXEDBLOCKSIZE
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  sizeX=blkLimitsGC(HIGH,IAXIS)
  sizeY=blkLimitsGC(HIGH,JAXIS)
  sizeZ=blkLimitsGC(HIGH,KAXIS)
  allocate(xCenter(sizeX))
  allocate(xLeft(sizeX))
  allocate(xRight(sizeX))
  allocate(yCenter(sizeY))
  allocate(yLeft(sizeY))
  allocate(yRight(sizeY))
  allocate(zCenter(sizeZ))
  allocate(zLeft(sizeZ))
  allocate(zRight(sizeZ))
#else
  sizeX=GRID_IHI_GC
  sizeY=GRID_JHI_GC
  sizeZ=GRID_KHI_GC
#endif
  zLeft = 0.
  zCenter = 0.
  zRight = 0.
  yLeft = 0.
  yCenter = 0.
  yRight = 0.
  xLeft = 0.
  xCenter = 0.
  xRight = 0.

  if (NDIM == 3) then
     call Grid_getCellCoords(KAXIS, blockId, CENTER, gcell,zCenter, sizeZ)
     call Grid_getCellCoords(KAXIS, blockId, LEFT_EDGE, gcell, zLeft, sizeZ)
     call Grid_getCellCoords(KAXIS, blockId, RIGHT_EDGE, gcell,zRight, sizeZ)
  end if
  if (NDIM > 1) then
     call Grid_getCellCoords(JAXIS, blockId, CENTER, gcell,yCenter, sizeY)
     call Grid_getCellCoords(JAXIS, blockId, LEFT_EDGE, gcell, yLeft, sizeY)
     call Grid_getCellCoords(JAXIS, blockId, RIGHT_EDGE, gcell,yRight, sizeY)
  end if

  call Grid_getCellCoords(IAXIS, blockId, CENTER, gcell, xCenter, sizeX)
  call Grid_getCellCoords(IAXIS, blockId, LEFT_EDGE, gcell, xLeft, sizeX)
  call Grid_getCellCoords(IAXIS, blockId, RIGHT_EDGE, gcell, xRight, sizeX)

  j = pos(1); k=pos(2)
  grav(1:numCells) = 0.
  do i = 1, numCells
        
     select case (sweepDir)

     case (SWEEP_X)
        call grv_accelOneZone(xCenter (i), yCenter(j), zCenter(k), gc)
        call grv_accelOneZone(xLeft(i), yCenter(j), zCenter(k), gl)
        call grv_accelOneZone(xRight(i), yCenter(j), zCenter(k), gr)

     case (SWEEP_Y)
        call grv_accelOneZone(xCenter(j), yCenter (i), zCenter(k), gc)
        call grv_accelOneZone(xCenter(j), yLeft(i), zCenter(k), gl)
        call grv_accelOneZone(xCenter(j), yRight(i), zCenter(k), gr)

     case (SWEEP_Z)
        call grv_accelOneZone(xCenter(j), yCenter(k), zCenter (i), gc)
        call grv_accelOneZone(xCenter(j), yCenter(k), zLeft(i), gl)
        call grv_accelOneZone(xCenter(j), yCenter(k), zRight(i), gr)
        
     end select

     grav_zone = (gl + 4.*gc + gr)/6.
     grav(i) = grav_zone(sweepDir)

  enddo
#ifndef FIXEDBLOCKSIZE
  deallocate(xCenter)
  deallocate(xLeft)
  deallocate(xRight)
  deallocate(yCenter)
  deallocate(yLeft)
  deallocate(yRight)
  deallocate(zCenter)
  deallocate(zLeft)
  deallocate(zRight)
#endif
!
!==============================================================================
!
  return
end subroutine Gravity_accelOneRow
