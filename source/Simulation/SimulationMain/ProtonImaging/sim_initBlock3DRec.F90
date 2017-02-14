!!****if* source/Simulation/SimulationMain/ProtonImaging/sim_initBlock3DRec
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
!!  Initializes the data (electric fields only -> 0) for a specified block needed
!!  to run the proton imaging. Specific routine for 3D rectangular geometries.
!!
!! ARGUMENTS
!!
!!  blockID : The block ID number to be initialized
!!
!!***

subroutine sim_initBlock3DRec (blockID)

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

  integer, dimension (1:2) :: dataSize
  integer, dimension (1:3) :: startPosition

  integer, dimension (LOW:HIGH,3) :: blkLimits
  integer, dimension (LOW:HIGH,3) :: blkLimitsGC

  real, allocatable :: cellEdgesX    (:)
  real, allocatable :: cellEdgesZ    (:)
  real, allocatable :: dataBlockEleX (:,:)
  real, allocatable :: dataBlockEleY (:,:)
  real, allocatable :: dataBlockEleZ (:,:)
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
  allocate (dataBlockEleX (1:nCellsX,1:nCellsZ))
  allocate (dataBlockEleY (1:nCellsX,1:nCellsZ))
  allocate (dataBlockEleZ (1:nCellsX,1:nCellsZ))

  call Grid_getCellCoords (IAXIS, blockID, FACES, includeGuardCells, cellEdgesX, nCellsX+1)
  call Grid_getCellCoords (KAXIS, blockID, FACES, includeGuardCells, cellEdgesZ, nCellsZ+1)
!
!
!    ...Set all electric components to zero.
!
!
  do i = 1,nCellsX
     do k = 1,nCellsZ
        dataBlockEleX (i,k) = 0.0
        dataBlockEleY (i,k) = 0.0
        dataBlockEleZ (i,k) = 0.0
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

  end do
!
!
!    ...Deallocate the intermediate arrays.
!
!
  deallocate (cellEdgesX   )
  deallocate (cellEdgesZ   )
  deallocate (dataBlockEleX)
  deallocate (dataBlockEleY)
  deallocate (dataBlockEleZ)
!
!
!    ...Ready!
!
!
  return
end subroutine sim_initBlock3DRec
