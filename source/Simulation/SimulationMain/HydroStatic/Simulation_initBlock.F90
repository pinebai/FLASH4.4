!!****if* source/Simulation/SimulationMain/HydroStatic/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!! 
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer :: blockId)
!!                       
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.  This version sets up the Hydrostatic boundary 
!!  condition test problem.
!!
!!
!! ARGUMENTS
!!
!!  blockId -        The number of the block to initialize
!!
!!
!!***

subroutine Simulation_initBlock(blockID)
  ! get the needed unit scope data
  use Grid_interface, ONLY : Grid_getBlkIndexLimits,Grid_getCellCoords,&
                              Grid_putPointData
  use Simulation_data, ONLY : sim_gamma, sim_xyzRef, sim_presRef, &
       sim_gravVector, sim_gravDirec, sim_gravConst, &
       sim_xyzRef, sim_densRef, sim_tempRef, &
       sim_molarMass, sim_gasconstant
  implicit none


! get all the constants
#include "constants.h"
#include "Flash.h"


! define arguments and indicate whether they are input or output
  integer, intent(in) :: blockID
  


! declare all local variables.
  integer :: i, j, k, n
  real :: xx, yy,  zz, gh

  ! arrays to hold coordinate information for the block
  real,allocatable, dimension(:) ::xCoord,yCoord,zCoord


  ! array to get integer indices defining the beginning and the end
  ! of a block.
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC


! the number of grid points along each dimension
  integer :: sizeX,sizeY,sizeZ

  integer, dimension(MDIM) :: axis
  integer,parameter :: dataSize = 1  ! for Grid put data function
  logical,parameter :: gcell = .true.


  ! these variables store the calculated initial values of physical
  ! variables a grid point at a time.
  real :: rhoZone, velxZone, velyZone, velzZone, presZone, &
       enerZone, ekinZone, tempZone

!  print*,' PE', myPE,' initializing block no.', blockID

  ! get the integer endpoints of the block in all dimensions
  ! the array blkLimits returns the interior end points
  ! whereas array blkLimitsGC returns endpoints including guardcells
  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)


! get the size along each dimension for allocation and then allocate
  sizeX = blkLimitsGC(HIGH,IAXIS)
  sizeY = blkLimitsGC(HIGH,JAXIS)
  sizeZ = blkLimitsGC(HIGH,KAXIS)
  allocate(xCoord(sizeX))
  allocate(yCoord(sizeY))
  allocate(zCoord(sizeZ))
  call Grid_getCellCoords(IAXIS, blockID, CENTER, gcell, xCoord, sizeX)
  call Grid_getCellCoords(JAXIS, blockID, CENTER, gcell, yCoord, sizeY)
  call Grid_getCellCoords(KAXIS, blockID, CENTER, gcell, zCoord, sizeZ)


  !-----------------------------------------------------------------------------
  ! loop over all of the zones in the current block and set the variables.
  !-----------------------------------------------------------------------------
  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     zz = zCoord(k) ! coordinates of the cell center in the z-direction

     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        yy = yCoord(j) ! center coordinates in the y-direction

        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           xx = xCoord(i)


           axis(IAXIS) = i   ! Get the position of the cell in the block
           axis(JAXIS) = j
           axis(KAXIS) = k

           ! Compute the gas energy and set the gamma-values
           ! needed for the equation of  state.

           tempZone = sim_tempRef

           gh = sim_gravVector(1) * (xx-sim_xyzRef) + &
                 sim_gravVector(2) * (yy-sim_xyzRef) + &
                 sim_gravVector(3) * (zz-sim_xyzRef)

!           print*,i,j,exp( sim_molarMass * gh / (sim_gasconstant * tempZone) )
           rhoZone = sim_densRef * exp( sim_molarMass * gh / (sim_gasconstant * tempZone) )
           presZone = sim_presRef * rhoZone / sim_densRef

           ekinZone = 0.0

           enerZone = presZone / (sim_gamma-1.)
           enerZone = enerZone / rhoZone
           enerZone = enerZone + ekinZone

           ! store the variables in the current zone via the Grid_putPointData method

           call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rhoZone)
           call Grid_putPointData(blockId, CENTER, PRES_VAR, EXTERIOR, axis, presZone)
           call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, axis, velxZone)
           call Grid_putPointData(blockId, CENTER, VELY_VAR, EXTERIOR, axis, velyZone)
           call Grid_putPointData(blockId, CENTER, VELZ_VAR, EXTERIOR, axis, velzZone)
           call Grid_putPointData(blockId, CENTER, ENER_VAR, EXTERIOR, axis, enerZone)
           call Grid_putPointData(blockId, CENTER, EINT_VAR, EXTERIOR, axis, enerZone)
           call Grid_putPointData(blockId, CENTER, TEMP_VAR, EXTERIOR, axis, tempZone)
           call Grid_putPointData(blockId, CENTER, GAME_VAR, EXTERIOR, axis, sim_gamma)
           call Grid_putPointData(blockId, CENTER, GAMC_VAR, EXTERIOR, axis, sim_gamma)

        end do
     end do
  end do

  deallocate(xCoord)
  deallocate(yCoord)
  deallocate(zCoord)


end subroutine Simulation_initBlock
