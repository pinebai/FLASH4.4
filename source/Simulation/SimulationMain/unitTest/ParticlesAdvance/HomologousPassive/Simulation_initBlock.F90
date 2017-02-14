!!****if* source/Simulation/SimulationMain/unitTest/ParticlesAdvance/HomologousPassive/Simulation_initBlock
!!
!! NAME
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!  Simulation_initBlock(  integer(in) :: blockID,
!!                         
!!
!! DESCRIPTION   
!!     Initialize fluid properties (density, pressure, velocity, etc.) 
!!          in a single block for the unitTest/Particles
!!     Set up uniform properties with constant density and pressure
!!     For velocities, y-vel and z-vel are zero.  In the x-direction,
!!       a constant value is initialized.  Half of the y domain can have a
!!       constant multiple of the x-velocity in the rest of the domain
!!
!! ARGUMENTS
!!      blockID:     integer(in)      the current block number to be filled
!!      
!!
!! PARAMETERS
!!
!!   sim_rho_amb    Gas Density:  Entire domain receives this ambient parameter
!!   sim_p_amb      Gas Pressure: Entire domain receives this ambient parameter
!!   sim_vx_amb     Gas x-velocity:  Dominant flow velocity throughout domain 
!!   sim_vx_multiplier   Half of the domain in y has x-velocity multiplied by this value
!!
!!***


subroutine Simulation_initBlock(blockID)

!============================================================================

 
  use Simulation_data, ONLY: sim_vx_amb, sim_rho_amb, sim_p_amb, sim_vx_multiplier
  use Logfile_interface, ONLY : Logfile_stamp
  use Driver_interface, ONLY : Driver_getSimTime
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_putPointData, Grid_getCellCoords
  use sim_interface
  implicit none 

#include "constants.h"
#include "Flash.h"


 integer, intent(IN) :: blockID

  integer               :: b, i, j, k, n, status 
  integer    :: ii, jj, kk
  integer    :: jhalf
  real       :: vx  ! velocity in x direction, may be increased in half of region with sim_vx_multiplier

  real,allocatable, dimension(:) ::xCenter,yCoord,zCoord

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer :: sizeX,sizeY,sizeZ
  integer, dimension(MDIM) :: axis

  real :: velxZone, velyZone, velzZone
  real :: tInit
  
  logical :: gcell = .true.

  ! Write a message to stdout describing the problem setup.
  call Logfile_stamp(                              & 
          "initializing for unitTest for Particles module", '[Simulation_initBlock]')
  write (*,*) "flash:  initializing for unitTest for Particles module"


  call Driver_getSimTime(tInit)
  print*,'tInit =',tInit
  
!----------------------------------------------------------------------------

  ! get cell coordinates for the block

! get the number of cells in each direction
!  blkLimitsGC includes guard cells, blkLimits ignores them.
  call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)

!----------------------------------------------------------------------------

! Initialize the flowfield.
!  velocity in x direction is vx_amb, velocity in other directions is zero
 
  ! Determine half of the y domain for multiplier
  jhalf = (blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS))/2 + blkLimitsGC(LOW,JAXIS)

  sizeX = blkLimitsGC(HIGH,IAXIS)
  sizeY = blkLimitsGC(HIGH,JAXIS)
  sizeZ = blkLimitsGC(HIGH,KAXIS)
  allocate(xCenter(sizeX))
  allocate(yCoord(sizeY))
  allocate(zCoord(sizeZ))
  xCenter = 0.0
  yCoord = 0.0
  zCoord = 0.0

  if (NDIM == 3) call Grid_getCellCoords&
                      (KAXIS, blockId, CENTER,gcell, zCoord, sizeZ)
  if (NDIM >= 2) call Grid_getCellCoords&
                      (JAXIS, blockId, CENTER,gcell, yCoord, sizeY)

  call Grid_getCellCoords(IAXIS, blockId, CENTER, gcell, xCenter, sizeX)

  ! Loop over (interior) cells in the block and initialize the variables.
  do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     axis(KAXIS) = k

     do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        axis(JAXIS) = j

        ! new option for changing the initial velocity with variable y
        vx = sim_vx_amb
        if (j .ge. jhalf) vx = sim_vx_amb*sim_vx_multiplier
        

        do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
           axis(IAXIS) = i

           !-------------------------------------------------------------------
           ! finally, fill the solution array
           !-------------------------------------------------------------------

           call Grid_putPointData(blockID, CENTER, DENS_VAR, EXTERIOR, axis, sim_rho_amb)
           call Grid_putPointData(blockID, CENTER, PRES_VAR, EXTERIOR, axis, sim_p_amb)
           call sim_grAnaGetVelComponent(velxZone,VELX_VAR,xCenter(i),yCoord(j),zCoord(k),tInit)
           call sim_grAnaGetVelComponent(velyZone,VELY_VAR,xCenter(i),yCoord(j),zCoord(k),tInit)
           call sim_grAnaGetVelComponent(velzZone,VELZ_VAR,xCenter(i),yCoord(j),zCoord(k),tInit)
           call Grid_putPointData(blockID, CENTER, VELX_VAR, EXTERIOR, axis, velxZone)
           call Grid_putPointData(blockID, CENTER, VELY_VAR, EXTERIOR, axis, velyZone)
           call Grid_putPointData(blockID, CENTER, VELZ_VAR, EXTERIOR, axis, velzZone)


        enddo  ! sweep over x
     enddo     ! sweep over y
  enddo        ! sweep over z

  !============================================================================

  print *,'Done with Simulation_initBlock'

  return
end subroutine Simulation_initBlock
