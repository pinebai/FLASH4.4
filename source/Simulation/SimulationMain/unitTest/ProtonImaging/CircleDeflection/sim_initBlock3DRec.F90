!!****if* source/Simulation/SimulationMain/unitTest/ProtonImaging/CircleDeflection/sim_initBlock3DRec
!!
!! NAME
!!
!!  sim_initBlock3DRec
!!
!! SYNOPSIS
!!
!!  sim_initBlock3DRec (integer (in) :: blockID)
!!
!! DESCRIPTION
!!
!!  Initializes the data (electric, magnetic fields and boundary indicator) for a
!!  specified block needed to run the proton imaging unit test. Specific routine
!!  for 3D rectangular geometries.
!!
!! ARGUMENTS
!!
!!  blockID : The block ID number to be initialized
!!
!!***

subroutine sim_initBlock3DRec (blockID)

  use Simulation_data,         ONLY : sim_clockwiseB,          &
                                      sim_magneticFluxDensity, &
                                      sim_xCenter,             &
                                      sim_zCenter

  use Driver_interface,        ONLY : Driver_abortFlash

  use Grid_interface,          ONLY : Grid_getBlkIndexLimits, &
                                      Grid_getCellCoords,     &
                                      Grid_putPlaneData

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent (in) :: blockID

  logical, save :: includeGuardCells = .false.

  integer  :: i,j,k
  integer  :: imin,imax
  integer  :: jmin,jmax
  integer  :: kmin,kmax
  integer  :: nCellsX, nCellsY, nCellsZ

  real     :: B,Bx,By,Bz
  real     :: cx,cz
  real     :: r,rx,rz
  real     :: s
  real     :: x1, x2, z1, z2

  integer, dimension (1:2) :: dataSize
  integer, dimension (1:3) :: startPosition

  integer, dimension (LOW:HIGH,3) :: blkLimits
  integer, dimension (LOW:HIGH,3) :: blkLimitsGC

  real, allocatable :: cellEdgesX    (:)
  real, allocatable :: cellEdgesZ    (:)
  real, allocatable :: dataBlockMagX (:,:)
  real, allocatable :: dataBlockMagY (:,:)
  real, allocatable :: dataBlockMagZ (:,:)
  real, allocatable :: dataBlockEleX (:,:)
  real, allocatable :: dataBlockEleY (:,:)
  real, allocatable :: dataBlockEleZ (:,:)
  real, allocatable :: dataBlockBdry (:,:)
!
!
!    ...Loop over all cells in current block and initialize the cells with
!       the needed data.
!
!
  call Grid_getBlkIndexLimits (blockID,    &
                               blkLimits,  &
                               blkLimitsGC )

  imin = blkLimits (LOW ,IAXIS)
  imax = blkLimits (HIGH,IAXIS)
  jmin = blkLimits (LOW ,JAXIS)
  jmax = blkLimits (HIGH,JAXIS)
  kmin = blkLimits (LOW ,KAXIS)
  kmax = blkLimits (HIGH,KAXIS)

  nCellsX = imax - imin + 1
  nCellsY = jmax - jmin + 1
  nCellsZ = kmax - kmin + 1

  allocate (cellEdgesX    (nCellsX+1))
  allocate (cellEdgesZ    (nCellsZ+1))
  allocate (dataBlockMagX (1:nCellsX,1:nCellsZ))
  allocate (dataBlockMagY (1:nCellsX,1:nCellsZ))
  allocate (dataBlockMagZ (1:nCellsX,1:nCellsZ))
  allocate (dataBlockEleX (1:nCellsX,1:nCellsZ))
  allocate (dataBlockEleY (1:nCellsX,1:nCellsZ))
  allocate (dataBlockEleZ (1:nCellsX,1:nCellsZ))
  allocate (dataBlockBdry (1:nCellsX,1:nCellsZ))

  call Grid_getCellCoords (IAXIS, blockID, FACES, includeGuardCells, cellEdgesX, nCellsX+1)
  call Grid_getCellCoords (KAXIS, blockID, FACES, includeGuardCells, cellEdgesZ, nCellsZ+1)
!
!
!    ...Currently sets all electric components to zero and the boundary indicator
!       to non-boundary (-1.0).
!
!
  do i = 1,nCellsX
     do k = 1,nCellsZ
        dataBlockEleX (i,k) = 0.0
        dataBlockEleY (i,k) = 0.0
        dataBlockEleZ (i,k) = 0.0
        dataBlockBdry (i,k) = -1.0
     end do
  end do
!
!
!    ...Set the magnetic components, such that the B vector is oriented tangentially at each
!       radial point. The orientation can be clockwise or anticlockwise along the radial
!       center.
!
!
  cx = sim_xCenter
  cz = sim_zCenter
  B  = sim_magneticFluxDensity / sqrt (4.0 * PI)  ! to conform with FLASH B's

  if (sim_clockwiseB) then
      s = 1.0
  else
      s = -1.0
  end if

  do i = 1,nCellsX

     x1 = cellEdgesX (i)
     x2 = cellEdgesX (i+1)
     rx  = 0.5 * (x1 + x2) - cx     ! local radial x-component from center

     do k = 1,nCellsZ

        z1 = cellEdgesZ (k)
        z2 = cellEdgesZ (k+1)
        rz = 0.5 * (z1 + z2) - cz   ! local radial z-component from center

        r = sqrt (rx * rx + rz * rz)

!        if ((rx == 0.0 .and. rz == 0.0) .or. (r < 0.1) .or. (r > 0.15)) then
        if (rx == 0.0 .and. rz == 0.0) then
            Bx = 0.0
            By = 0.0
            Bz = 0.0
        else if (abs (rx) <= abs (rz)) then
            r  = rx / rz
            Bx = s * sign (1.0,rz) * B  / sqrt (1.0 + r * r)
            Bz = - r * Bx
            By = 0.0
        else
            r  = rz / rx
            Bz = - s * sign (1.0,rx) * B / sqrt (1.0 + r * r)
            Bx = - r * Bz
            By = 0.0
        end if

        dataBlockMagX (i,k) = Bx
        dataBlockMagY (i,k) = By
        dataBlockMagZ (i,k) = Bz

     end do
  end do
!
!
!    ...Put the individual xz data planes.
!
!
  dataSize (1) = nCellsX
  dataSize (2) = nCellsZ

  do j = 1,nCellsY

     startPosition (1) = 1
     startPosition (2) = j
     startPosition (3) = 1

     call Grid_putPlaneData (blockID,       &
                             CENTER,        &
                             MAGX_VAR,      &
                             INTERIOR,      &
                             XZPLANE,       &
                             startPosition, &
                             dataBlockMagX, &
                             dataSize       )

     call Grid_putPlaneData (blockID,       &
                             CENTER,        &
                             MAGY_VAR,      &
                             INTERIOR,      &
                             XZPLANE,       &
                             startPosition, &
                             dataBlockMagY, &
                             dataSize       )

     call Grid_putPlaneData (blockID,       &
                             CENTER,        &
                             MAGZ_VAR,      &
                             INTERIOR,      &
                             XZPLANE,       &
                             startPosition, &
                             dataBlockMagZ, &
                             dataSize       )

     call Grid_putPlaneData (blockID,       &
                             CENTER,        &
                             ELEX_VAR,      &
                             INTERIOR,      &
                             XZPLANE,       &
                             startPosition, &
                             dataBlockEleX, &
                             dataSize       )

     call Grid_putPlaneData (blockID,       &
                             CENTER,        &
                             ELEY_VAR,      &
                             INTERIOR,      &
                             XZPLANE,       &
                             startPosition, &
                             dataBlockEleY, &
                             dataSize       )

     call Grid_putPlaneData (blockID,       &
                             CENTER,        &
                             ELEZ_VAR,      &
                             INTERIOR,      &
                             XZPLANE,       &
                             startPosition, &
                             dataBlockEleZ, &
                             dataSize       )

     call Grid_putPlaneData (blockID,       &
                             CENTER,        &
                             BDRY_VAR,      &
                             INTERIOR,      &
                             XZPLANE,       &
                             startPosition, &
                             dataBlockBdry, &
                             dataSize       )
  end do
!
!
!    ...Deallocate the intermediate arrays.
!
!
  deallocate (cellEdgesX   )
  deallocate (cellEdgesZ   )
  deallocate (dataBlockMagX)
  deallocate (dataBlockMagY)
  deallocate (dataBlockMagZ)
  deallocate (dataBlockEleX)
  deallocate (dataBlockEleY)
  deallocate (dataBlockEleZ)
  deallocate (dataBlockBdry)
!
!
!    ...Ready!
!
!
  return
end subroutine sim_initBlock3DRec
