!!****if* source/Simulation/SimulationMain/unitTest/Laser_quadraticTube/testII/sim_initBlock3DRec
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
!!  Initializes the data (density and temperatures) for a specified block needed to run
!!  the the laser quadratic tube unit test. Specific routine for 3D rectangular geometries.
!!
!! ARGUMENTS
!!
!!  blockID : The block ID number to be initialized
!!
!!***

subroutine sim_initBlock3DRec (blockID)

  use EnergyDeposition_data,   ONLY : ed_Avogadro

  use Simulation_data,         ONLY : sim_A,  &
                                      sim_nc, &
                                      sim_nw, &
                                      sim_Tw, &
                                      sim_xc, &
                                      sim_xw, &
                                      sim_zc, &
                                      sim_zw

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

  real     :: Athird
  real     :: cellDensity
  real     :: cellTele
  real     :: cellTemp
  real     :: dx, dz
  real     :: ne, Te
  real     :: third, twothirds
  real     :: x1, x2, xm
  real     :: z1, z2, zm

  integer, dimension (1:2) :: dataSize
  integer, dimension (1:3) :: startPosition

  integer, dimension (LOW:HIGH,3) :: blkLimits
  integer, dimension (LOW:HIGH,3) :: blkLimitsGC

  real, allocatable :: cellEdgesX    (:)
  real, allocatable :: cellEdgesZ    (:)
  real, allocatable :: dataBlockDens (:,:)
  real, allocatable :: dataBlockTele (:,:)
  real, allocatable :: dataBlockTemp (:,:)
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
  allocate (dataBlockDens (1:nCellsX,1:nCellsZ))
  allocate (dataBlockTele (1:nCellsX,1:nCellsZ))
  allocate (dataBlockTemp (1:nCellsX,1:nCellsZ))

  call Grid_getCellCoords (IAXIS, blockID, FACES, includeGuardCells, cellEdgesX, nCellsX+1)
  call Grid_getCellCoords (KAXIS, blockID, FACES, includeGuardCells, cellEdgesZ, nCellsZ+1)
!
!
!    ...Fill in the cell mass density in each cell, such that a quadratic circular tube
!       electron density is obtained by the laser energy deposition code. The quadratic
!       circular tube is defined for the 3D cartesian geometry in the following way:
!
!           i) base of the tube   --> the domain xz-plane (square)
!          ii) length of the tube --> length of domain y-axis
!
!       The goal is to provide a density for each cell, such that the computed electron
!       number density will correspond to the quadratic tube:
!
!                 ne (x,z) = nw + A[(x - xw)^2 + (z - zw)^2]     A = (nc - nw) / (xc - xw)^2
!
!       To set the 'ne' for each cell, we integrate over all ne (x,z) inside the cell
!       square and divide by the cell xz-area.
!
!
!                       z-axis |
!                              |
!                              |
!                             z2        -----------
!                              |       |           |
!                              |       |     xz    |
!                              |       |   plane   |
!                              |       |    cell   |
!                              |       |           |
!                             z1        -----------
!                              |
!                              |
!                              |
!                               ------x1-----------x2-----------> x-axis
!
!
!       The result is:
!
!            ne (cell) = nw + A[(xm - xw)^2 + (zm - zw)^2] + (A/3)[dx^2 + dz^2]
!
!       where:
!
!                              xm = (x1 + x2) / 2      (cell midpoint)
!                              zm = (z1 + z2) / 2      (cell midpoint)
!                              dx = (x1 - x2) / 2
!                              dz = (z1 - z2) / 2
!
!       Inside the energy deposition code, the 'ne' are calculated as follows:
!
!                ne (cell) = cellDensity (cell) * nAvogadro * Zbar / Abar
!
!       The values of Zbar and Abar are set equal to 1. for all cells, hence we are left with:
!
!                ne (cell) = cellDensity (cell) * nAvogadro
!
!       The required equation for each cell density becomes then:
!
!                cellDensity (cell) = (1 / nAvogadro) * ne (cell)
!
!       These are the densities that will be stored. The electron temperature is obtained as:
!
!                   Te (cell) = [ne (cell) / nw]^(2/3)
!
!
!
  third     = 1.0 / 3.0
  twothirds = third + third
  Athird    = sim_A * third

  do i = 1,nCellsX

     x1 = cellEdgesX (i)
     x2 = cellEdgesX (i+1)

     xm = 0.5 * (x1 + x2)
     dx = 0.5 * (x1 - x2)

     do k = 1,nCellsZ

        z1 = cellEdgesZ (k)
        z2 = cellEdgesZ (k+1)

        zm = 0.5 * (z1 + z2)
        dz = 0.5 * (z1 - z2)

        ne = sim_nw + sim_A * ((xm - sim_xw) ** 2 + (zm - sim_zw) ** 2) + Athird * (dx * dx + dz * dz)
        Te = sim_Tw * ((ne / sim_nw) ** twothirds)

        cellDensity = ne / ed_Avogadro
        cellTele    = Te
        cellTemp    = cellTele

        dataBlockDens (i,k) = cellDensity
        dataBlockTele (i,k) = cellTele
        dataBlockTemp (i,k) = cellTemp

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
                             DENS_VAR,      &
                             INTERIOR,      &
                             XZPLANE,       &
                             startPosition, &
                             dataBlockDens, &
                             dataSize       )

     call Grid_putPlaneData (blockID,       &
                             CENTER,        &
                             TELE_VAR,      &
                             INTERIOR,      &
                             XZPLANE,       &
                             startPosition, &
                             dataBlockTele, &
                             dataSize       )

     call Grid_putPlaneData (blockID,       &
                             CENTER,        &
                             TEMP_VAR,      &
                             INTERIOR,      &
                             XZPLANE,       &
                             startPosition, &
                             dataBlockTemp, &
                             dataSize       )

  end do
!
!
!    ...Deallocate the intermediate arrays.
!
!
  deallocate (cellEdgesX   )
  deallocate (cellEdgesZ   )
  deallocate (dataBlockDens)
  deallocate (dataBlockTele)
  deallocate (dataBlockTemp)
!
!
!    ...Ready!
!
!
  return
end subroutine sim_initBlock3DRec
