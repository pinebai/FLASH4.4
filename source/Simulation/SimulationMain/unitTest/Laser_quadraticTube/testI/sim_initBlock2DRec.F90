!!****if* source/Simulation/SimulationMain/unitTest/Laser_quadraticTube/testI/sim_initBlock2DRec
!!
!! NAME
!!
!!  sim_initBlock2DRec
!!
!! SYNOPSIS
!!
!!  sim_initBlock2DRec (integer (in) :: blockID)
!!
!! DESCRIPTION
!!
!!  Initializes the data (density and temperatures) for a specified block needed to run
!!  the the laser quadratic tube unit test. Specific routine for 2D rectangular geometries.
!!
!! ARGUMENTS
!!
!!  blockID : The block ID number to be initialized
!!
!!***

subroutine sim_initBlock2DRec (blockID)

  use EnergyDeposition_data,   ONLY : ed_Avogadro

  use Simulation_data,         ONLY : sim_A,                &
                                      sim_nc,               &
                                      sim_nw,               &
                                      sim_Tw,               &
                                      sim_xc,               &
                                      sim_yc,               &
                                      sim_xw,               &
                                      sim_yw,               &
                                      sim_lasersOrientation

  use Driver_interface,        ONLY : Driver_abortFlash

  use Grid_interface,          ONLY : Grid_getBlkIndexLimits, &
                                      Grid_getCellCoords,     &
                                      Grid_putRowData

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent (in) :: blockID

  logical, save :: includeGuardCells = .false.

  integer  :: dataSize
  integer  :: i,j
  integer  :: imin,imax
  integer  :: jmin,jmax
  integer  :: nCellsX, nCellsY

  real     :: Athird
  real     :: cellDensity
  real     :: cellTele
  real     :: cellTemp
  real     :: dx, dy
  real     :: ne, Te
  real     :: third, twothirds
  real     :: x1, x2, xm
  real     :: y1, y2, ym

  integer, dimension (1:3) :: startPosition

  integer, dimension (LOW:HIGH,3) :: blkLimits
  integer, dimension (LOW:HIGH,3) :: blkLimitsGC

  real, allocatable :: cellEdges     (:)
  real, allocatable :: dataBlockDens (:)
  real, allocatable :: dataBlockTele (:)
  real, allocatable :: dataBlockTemp (:)
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

  nCellsX = imax - imin + 1
  nCellsY = jmax - jmin + 1

  third     = 1.0 / 3.0
  twothirds = third + third
  Athird    = sim_A * third
!
!
!    ...Select according to lasers orientation. Shown are details of computation for lasers
!       orientation on to the X domain side.
!
!       Fill in the cell mass density in each cell, such that a quadratic 2D tube
!       electron density is obtained by the laser energy deposition code. The quadratic
!       2D tube is defined for the 2D cartesian geometry in the following way:
!
!           i) base of the tube   --> the domain x-line
!          ii) length of the tube --> length of domain y-axis
!
!       The goal is to provide a density for each cell, such that the computed electron
!       number density will correspond to the quadratic 2D tube:
!
!                 ne (x) = nw + A(x - xw)^2     A = (nc - nw) / (xc - xw)^2
!
!       To set the 'ne' for each cell, we integrate over all ne (x) within the cell's line
!       and divide by the length of the line. The result is:
!
!            ne (cell) = nw + A(xm - xw)^2 + (A/3)dx^2
!
!       where:
!
!                              xm = (x1 + x2) / 2      (cell midpoint)
!                              dx = (x1 - x2) / 2
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
  select case (sim_lasersOrientation)

    case ('X')

      allocate (cellEdges     (nCellsX+1))
      allocate (dataBlockDens (1:nCellsY))
      allocate (dataBlockTele (1:nCellsY))
      allocate (dataBlockTemp (1:nCellsY))

      call Grid_getCellCoords (IAXIS, blockID, FACES, includeGuardCells, cellEdges, nCellsX+1)

      do i = 1,nCellsX

         x1 = cellEdges (i)
         x2 = cellEdges (i+1)

         xm = 0.5 * (x1 + x2)
         dx = 0.5 * (x1 - x2)

         ne = sim_nw + sim_A * (xm - sim_xw) ** 2 + Athird * dx * dx
         Te = sim_Tw * ((ne / sim_nw) ** twothirds)

         cellDensity = ne / ed_Avogadro
         cellTele    = Te
         cellTemp    = cellTele

         do j = 1,nCellsY                         ! prepare the individual data rows
            dataBlockDens (j) = cellDensity
            dataBlockTele (j) = cellTele
            dataBlockTemp (j) = cellTemp
         end do

         dataSize = nCellsY         ! transfer the data rows on to the grid

         startPosition (1) = i
         startPosition (2) = 1
         startPosition (3) = 1

         call Grid_putRowData (blockID,       &
                               CENTER,        &
                               DENS_VAR,      &
                               INTERIOR,      &
                               JAXIS,         &
                               startPosition, &
                               dataBlockDens, &
                               dataSize       )

         call Grid_putRowData (blockID,       &
                               CENTER,        &
                               TELE_VAR,      &
                               INTERIOR,      &
                               JAXIS,         &
                               startPosition, &
                               dataBlockTele, &
                               dataSize       )

         call Grid_putRowData (blockID,       &
                               CENTER,        &
                               TEMP_VAR,      &
                               INTERIOR,      &
                               JAXIS,         &
                               startPosition, &
                               dataBlockTemp, &
                               dataSize       )
      end do

      deallocate (cellEdges)
      deallocate (dataBlockDens)
      deallocate (dataBlockTele)
      deallocate (dataBlockTemp)

    case ('Y')

      allocate (cellEdges     (nCellsY+1))
      allocate (dataBlockDens (1:nCellsX))
      allocate (dataBlockTele (1:nCellsX))
      allocate (dataBlockTemp (1:nCellsX))

      call Grid_getCellCoords (JAXIS, blockID, FACES, includeGuardCells, cellEdges, nCellsY+1)

      do j = 1,nCellsY

         y1 = cellEdges (j)
         y2 = cellEdges (j+1)

         ym = 0.5 * (y1 + y2)
         dy = 0.5 * (y1 - y2)

         ne = sim_nw + sim_A * (ym - sim_yw) ** 2 + Athird * dy * dy
         Te = sim_Tw * ((ne / sim_nw) ** twothirds)

         cellDensity = ne / ed_Avogadro
         cellTele    = Te
         cellTemp    = cellTele

         do i = 1,nCellsX                         ! prepare the individual data rows
            dataBlockDens (i) = cellDensity
            dataBlockTele (i) = cellTele
            dataBlockTemp (i) = cellTemp
         end do

         dataSize = nCellsX         ! transfer the data rows on to the grid

         startPosition (1) = 1
         startPosition (2) = j
         startPosition (3) = 1

         call Grid_putRowData (blockID,       &
                               CENTER,        &
                               DENS_VAR,      &
                               INTERIOR,      &
                               IAXIS,         &
                               startPosition, &
                               dataBlockDens, &
                               dataSize       )

         call Grid_putRowData (blockID,       &
                               CENTER,        &
                               TELE_VAR,      &
                               INTERIOR,      &
                               IAXIS,         &
                               startPosition, &
                               dataBlockTele, &
                               dataSize       )

         call Grid_putRowData (blockID,       &
                               CENTER,        &
                               TEMP_VAR,      &
                               INTERIOR,      &
                               IAXIS,         &
                               startPosition, &
                               dataBlockTemp, &
                               dataSize       )
      end do

      deallocate (cellEdges)
      deallocate (dataBlockDens)
      deallocate (dataBlockTele)
      deallocate (dataBlockTemp)

    case default

      call Driver_abortFlash ('Laser quadratic tube unit test I: Invalid lasers orientation!')

  end select
!
!
!    ...Ready!
!
!
  return
end subroutine sim_initBlock2DRec
