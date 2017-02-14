!!****if* source/Simulation/SimulationMain/unitTest/Opacity/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock (integer (in) :: blockID)
!!
!! DESCRIPTION
!!
!!  Initializes the data (mass density, mass fraction and electron temperature) for
!!  a specified block needed to test the opacity unit.
!!
!! ARGUMENTS
!!
!!  blockID  :   The number of the block to be initialized
!!
!!***

subroutine Simulation_initBlock (blockID)

  use Driver_interface,     ONLY : Driver_abortFlash
  use Grid_interface,       ONLY : Grid_getBlkIndexLimits, &
                                   Grid_putRowData

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent(in) :: blockID

  integer  :: cell
  integer  :: mfPos
  integer  :: species
  integer  :: totalCells

  integer, dimension   (MDIM) :: startingPos
  integer, dimension (2,MDIM) :: blkLimits
  integer, dimension (2,MDIM) :: blkLimitsGC

  real :: minMassDensity
  real :: minMassFraction
  real :: minTemperature
  real :: massDensityStep
  real :: massFractionStep
  real :: temperatureStep

  real, dimension(:),   allocatable :: electronTemperature
  real, dimension(:),   allocatable :: massDensity
  real, dimension(:,:), allocatable :: massFraction
!
!
!    ...Safety net.
!
!
  if (NDIM /= 1) then
      call Driver_abortFlash ('[Simulation_initBlock] ERROR: # of dimensions > 1')
  end if
!
!
!    ...Loop over all cells in current block and initialize the cells with
!       the needed data. Initialize only the inner block cells, not the guard
!       cells. Work on a row-by-row basis.
!
!
  startingPos (1) = 1
  startingPos (2) = 1
  startingPos (3) = 1

  call Grid_getBlkIndexLimits (blockID,    &
                               blkLimits,  &
                               blkLimitsGC )

  totalCells = blkLimits (HIGH,IAXIS) - blkLimits (LOW,IAXIS) + 1

  allocate (massDensity         (totalCells))
  allocate (electronTemperature (totalCells))
  allocate (massFraction        (totalCells,NSPECIES))

  minMassDensity   = 1.E-6                     ! in g/cm^3
  minMassFraction  = 1./real(totalCells)       ! no units
  minTemperature   = 100.                      ! in K
  massDensityStep  = 10.
  massFractionStep = 1./real(totalCells)
  temperatureStep  = 4.
!
!
!    ...Set the mass density, mass fraction and electron temperature of each cell.
!       Currently each species mass fraction is set equal to 1/NSPECIES throughout
!       all cells.
!
!
  do cell = 1,totalCells
     massDensity         (cell) = minMassDensity  * massDensityStep ** (cell - 1)
!     massFraction        (cell) = minMassFraction + massFractionStep * (cell - 1)
     electronTemperature (cell) = minTemperature  * temperatureStep ** (cell - 1)

     do species = 1,NSPECIES
        massFraction (cell,species) = 1./real(NSPECIES)
     end do
  enddo
!
!
!    ...Store all the cell data.
!
!
  mfPos = SPECIES_BEGIN

  call Grid_putRowData (blockID,CENTER,DENS_VAR,INTERIOR,IAXIS,startingPos,massDensity,        totalCells)
  call Grid_putRowData (blockID,CENTER,TELE_VAR,INTERIOR,IAXIS,startingPos,electronTemperature,totalCells)

  do species = 1,NSPECIES
     mfPos = SPECIES_BEGIN + species - 1
     call Grid_putRowData (blockID,CENTER,mfPos,INTERIOR,IAXIS,startingPos,massFraction (:,species),totalCells)
  end do
!
!
!    ...Remove auxilliary arrays.
!
!
  deallocate (massDensity)
  deallocate (massFraction)
  deallocate (electronTemperature)
!
!
!    ...Ready!
!
!
  return
end subroutine Simulation_initBlock
