!!****if* source/Simulation/SimulationMain/unitTest/Laser_quadraticTube/RingExpPotential/sim_initBlock3DRec
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
!!  the laser ring exponential potential unit test. Specific routine for 3D rectangular
!!  geometries.
!!
!! ARGUMENTS
!!
!!  blockID : The block ID number to be initialized
!!
!!***

subroutine sim_initBlock3DRec (blockID)

  use EnergyDeposition_data,   ONLY : ed_Avogadro

  use Simulation_data,         ONLY : sim_nc,          &
                                      sim_Te0,         &
                                      sim_xTubeCenter, &
                                      sim_zTubeCenter, &
                                      sim_alpha

  use Driver_interface,        ONLY : Driver_abortFlash

  use Grid_interface,          ONLY : Grid_getBlkIndexLimits, &
                                      Grid_getCellCoords,     &
                                      Grid_putPlaneData

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent (in) :: blockID

  logical, save :: includeGuardCells = .false.

  integer  :: alpha
  integer  :: i,j,k
  integer  :: imin,imax
  integer  :: jmin,jmax
  integer  :: kmin,kmax
  integer  :: nCellsX, nCellsY, nCellsZ

  real     :: cellDensity
  real     :: cellTele
  real     :: ne, nc
  real     :: r, x, z
  real     :: Te, Te0
  real     :: x1, x2, z1, z2
  real     :: xc, zc

  integer, dimension (1:2) :: dataSize
  integer, dimension (1:3) :: startPosition

  integer, dimension (LOW:HIGH,3) :: blkLimits
  integer, dimension (LOW:HIGH,3) :: blkLimitsGC

  real, allocatable :: cellEdgesX    (:)
  real, allocatable :: cellEdgesZ    (:)
  real, allocatable :: dataBlockDens (:,:)
  real, allocatable :: dataBlockTele (:,:)
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

  call Grid_getCellCoords (IAXIS, blockID, FACES, includeGuardCells, cellEdgesX, nCellsX+1)
  call Grid_getCellCoords (KAXIS, blockID, FACES, includeGuardCells, cellEdgesZ, nCellsZ+1)
!
!
!    ...Fill in the cell mass density in each cell, such that a tube electron density is
!       obtained for a exponential potential by the laser energy deposition code. The circular
!       tube is defined for the 3D cartesian geometry in the following way:
!
!           i) base of the tube   --> the domain xz-plane (square)
!          ii) length of the tube --> length of domain y-axis
!
!       The goal is to provide a density for each cell, such that the computed electron
!       number density will correspond to:
!
!                             ne (r) = (nc / 4) * (1 + r ** alpha)
!
!       where
!
!                       nc = critical electron number density
!                        r = radius from center of tube
!                    alpha = potential exponent
!
!       The radius of a cell is calculated as
!
!                   r = sqrt ( (x - xc)^2 + (z - zc)^2 )
!
!       where (xc,zc) is the tube's center location and (x,z) is the location of the cell
!       center in the xz-plane:
!
!                              x = (x1 + x2) / 2      x1,x2 = cell x-axis edges
!                              z = (z1 + z2) / 2      z1,z2 = cell z-axis edges
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
!       These are the densities that will be stored. The electron temperature is evaluated
!       as (Te0 = electron temperature at tube center):
!
!                          Te (cell) = Te0 * (ne (cell) / ne0) ** (2/3)
!
!                                    = Te0 * (4 * ne (cell) / nc) ** (2/3)   (because ne0 = nc / 4)
!
!
!
  xc    = sim_xTubeCenter
  zc    = sim_zTubeCenter
  nc    = sim_nc
  Te0   = sim_Te0
  alpha = sim_alpha

  do i = 1,nCellsX

     x1 = cellEdgesX (i)
     x2 = cellEdgesX (i+1)
     x  = 0.5 * (x1 + x2)

     do k = 1,nCellsZ

        z1 = cellEdgesZ (k)
        z2 = cellEdgesZ (k+1)
        z  = 0.5 * (z1 + z2)

        r = sqrt ( (x - xc) ** 2 + (z - zc) ** 2 )

        ne = 0.25 * nc * (1.0 + r ** alpha)
        Te = Te0 * (4.0 * ne / nc) ** (2.0 / 3.0)

        cellDensity = ne / ed_Avogadro
        cellTele    = Te

        dataBlockDens (i,k) = cellDensity
        dataBlockTele (i,k) = cellTele

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
!
!
!    ...Ready!
!
!
  return
end subroutine sim_initBlock3DRec
