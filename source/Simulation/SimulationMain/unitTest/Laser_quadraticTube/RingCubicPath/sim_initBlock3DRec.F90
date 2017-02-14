!!****if* source/Simulation/SimulationMain/unitTest/Laser_quadraticTube/RingCubicPath/sim_initBlock3DRec
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
!!  the laser ring cubic path unit test. Specific routine for 3D rectangular geometries.
!!
!! ARGUMENTS
!!
!!  blockID : The block ID number to be initialized
!!
!!***

subroutine sim_initBlock3DRec (blockID)

  use EnergyDeposition_data,   ONLY : ed_Avogadro

  use Simulation_data,         ONLY : sim_nc,     &
                                      sim_nw,     &
                                      sim_Tw,     &
                                      sim_xc,     &
                                      sim_xw,     &
                                      sim_zc,     &
                                      sim_zw,     &
                                      sim_yfocal, &
                                      sim_R

  use Driver_interface,        ONLY : Driver_abortFlash

  use Grid_interface,          ONLY : Grid_getBlkIndexLimits, &
                                      Grid_getCellCoords,     &
                                      Grid_putPlaneData

  use ed_extractBeamsData,     ONLY : ed_extractBeamData

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

  real     :: cellDensity
  real     :: cellTele
  real     :: cellTemp
  real     :: cos1T3, cos2T3, sin1T3, sin2T3
  real     :: ne, Te
  real     :: nw, Tw
  real     :: r0, rFunc
  real     :: third, twothirds
  real     :: p, q, r, s, v, x, z
  real     :: posTerm, negTerm
  real     :: T
  real     :: x1, x2, z1, z2
  real     :: xw, zw
  real     :: yF

  real     :: rdiff,rmin


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
!    ...Fill in the cell mass density in each cell, such that a tube electron density is
!       obtained for a cubic ray path by the laser energy deposition code. The circular
!       tube is defined for the 3D cartesian geometry in the following way:
!
!           i) base of the tube   --> the domain xz-plane (square)
!          ii) length of the tube --> length of domain y-axis
!
!       The goal is to provide a density for each cell, such that the computed electron
!       number density will correspond to:
!
!                 ne (r) = nw * {1 + (3v/2yF)^2 * rFunc}    v = ray velocity in units of c
!                                                          yF = focal point along y axis
!
!       where:
!
!             rFunc =          2 * r0 * sqrt (r * (r0 - r)) * cos (T/3)
!                     +               (r0^2 + 2(r0 - 2r)^2) * cos (2T/3)
!                     +                  3 * r0 * (r0 - 2r) * sin (T/3)
!                     - 4 * (r0 - 2r) * sqrt (r * (r0 - r)) * sin (2T/3)
!                     - 3 * r0^2
!
!       and
!
!                 T = arctan ( (r0 - 2r) / 2 * sqrt (r * (r0 - r))
!
!       and
!
!                 r = sqrt ( (x - xw)^2 + (z - zw)^2 )
!
!       where (x,z) is the location of the cell center in the xZ-plane:
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
!       These are the densities that will be stored. The electron temperature is obtained as:
!
!                   Te (cell) = Tw * [ne (cell) / nw]^(2/3)
!
!
!
  third     = 1.0 / 3.0
  twothirds = third + third

  nw = sim_nw
  Tw = sim_Tw
  xw = sim_xw
  zw = sim_zw
  r0 = sim_R
  yF = 10.0
!  yF = sim_yfocal
!
!
!    ...Get the ray injection velocity (in units of c).
!
!
  call ed_extractBeamdata    (1,"initialRaySpeed",    v)

  s = (3 * v) / (2 * yF)
  s = s * s

  do i = 1,nCellsX

     x1 = cellEdgesX (i)
     x2 = cellEdgesX (i+1)
     x  = 0.5 * (x1 + x2)

     do k = 1,nCellsZ

        z1 = cellEdgesZ (k)
        z2 = cellEdgesZ (k+1)
        z  = 0.5 * (z1 + z2)

        r = sqrt ( (x - xw) ** 2 + (z - zw) ** 2 )

        if (r >= r0) then
            ne = nw * (1.0 + (12.0 * v * v) / (yF * yF) * r0 * (r - r0))
        else
            p = sqrt (r * (r0 - r))
            q = r0 - r - r
            T = atan (q / (p + p))

            cos1T3 = cos (T * third)
            cos2T3 = cos (T * twothirds)
            sin1T3 = sin (T * third)
            sin2T3 = sin (T * twothirds)

            posTerm = (r0 + r0) * p * cos1T3 + (r0 * r0 + 2 * q * q) * cos2T3 + 3 * r0 * q * sin1T3
            negTerm = 3 * r0 * r0 + 4 * p * q * sin2T3

            rFunc = posTerm - negTerm

            ne = nw * (1.0 + s * rFunc)
        end if

        Te = Tw * ((ne / nw) ** twothirds)

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

!  r = 1.50276437981807
!  p = sqrt (r * (r0 - r))
!  q = r0 - r - r
!  T = atan (q / (p + p))
!  cos1T3 = cos (T * third)
!  cos2T3 = cos (T * twothirds)
!  sin1T3 = sin (T * third)
!  sin2T3 = sin (T * twothirds)
!  posTerm = (r0 + r0) * p * cos1T3 + (r0 * r0 + 2 * q * q) * cos2T3 + 3 * r0 * q * sin1T3
!  negTerm = 3 * r0 * r0 + 4 * p * q * sin2T3
!  rFunc = posTerm - negTerm
!  ne = nw * (1.0 + s * rFunc)
!  write (*,*) ' ne = ',ne
!  write (*,*) ' ------------------------- '
!  rmin = r0

