!!****if* source/Simulation/SimulationMain/radflaHD/CriticalRadShock/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(IN) :: blockID) 
!!                       
!!
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.
!! 
!! ARGUMENTS
!!
!!  blockID -        the number of the block to initialize
!!  
!!
!!
!!***



subroutine Simulation_initBlock(blockId)
  use Simulation_data
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
       Grid_getCellCoords, Grid_putPointData, Grid_getPointData
  use Driver_interface, ONLY: Driver_abortFlash, Driver_getMype
  use RadTrans_interface, ONLY: RadTrans_mgdEFromT
  use Eos_interface, ONLY : Eos_wrapped

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"    

  integer, intent(in) :: blockId
  
  integer :: i, j, k
  integer :: istat
  integer :: blkLimits(2, MDIM)
  integer :: blkLimitsGC(2, MDIM)
  integer :: axis(MDIM)
  integer :: sizeX, sizeY, sizeZ
  real, allocatable :: xcent(:), ycent(:), zcent(:)
  real,allocatable,dimension(:) :: xCoord,yCoord,zCoord

  real :: rho   ! Density
  real :: tgas  ! Gas temperature 
  real :: trad  ! Radiation temperature
  real :: velx  ! 1D velocity
  real :: tradActual
  real :: xctr, yctr, zctr
  real :: xdist, ydist, zdist
  real :: dist,distxy
  real :: xrad, yrad, zrad, ekin
  
  real :: gasConst = 8.31447E+07
  real :: gamma    = 5.0/3.0
  real :: clight   = 2.99792E+10  
  real :: L_ref  = 1.0
  real :: radConst = 7.56577E-15

  real :: abs_opac    = 1.0E6
  real :: kappa       = 1.0  

  integer :: meshMe
  
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  
  sizeX = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  sizeY = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
  sizeZ = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

  allocate(xCoord(sizeX),stat=istat)
  allocate(yCoord(sizeY),stat=istat)
  allocate(zCoord(sizeZ),stat=istat)

  xCoord = 0.0
  yCoord = 0.0
  zCoord = 0.0

  call Grid_getCellCoords(KAXIS,blockID,CENTER,sim_gCell,zCoord,sizeZ)
  call Grid_getCellCoords(JAXIS,blockID,CENTER,sim_gCell,yCoord,sizeY)
  call Grid_getCellCoords(IAXIS,blockID,CENTER,sim_gCell,xCoord,sizeX)
  
  ! Calculate the coordinates of the center of the box
  !xctr = sim_xmin + (sim_xmax - sim_xmin)/2.0
  !yctr = sim_ymin + (sim_ymax - sim_ymin)/2.0
  !zctr = sim_zmin + (sim_zmax - sim_zmin)/2.0
 
  xrad = (sim_xmax - sim_xmin)/2.0

  ! Loop over cells and set the initial state
  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)          
           
           rho = sim_rho
           tgas = sim_tgas
           trad = sim_trad
           velx = sim_velx 
 !          ekin = 0.5*velx**2

           axis(IAXIS) = i
           axis(JAXIS) = j
           axis(KAXIS) = k
           
           call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rho)
           call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, axis, velx)
           call Grid_putPointData(blockId, CENTER, TEMP_VAR, EXTERIOR, axis, tgas)
           call Grid_putPointData(blockId, CENTER, TELE_VAR, EXTERIOR, axis, tgas)
           call Grid_putPointData(blockId, CENTER, TION_VAR, EXTERIOR, axis, 0.0)
#ifdef H1_SPEC
           call Grid_putPointData(blockId, CENTER, H1_SPEC, EXTERIOR, axis, 1.0)
#endif
#ifdef FLLM_VAR
           call Grid_putPointData(blockId, CENTER, FLLM_VAR, EXTERIOR, axis, 1.0)
#endif
           
           Call RadTrans_mgdEFromT(blockId, axis, trad, tradActual)
           call Grid_putPointData(blockId, CENTER, TRAD_VAR, EXTERIOR, axis, tradActual)
           
        enddo
     enddo
  enddo
    
  call Eos_wrapped(MODE_DENS_TEMP_GATHER, blkLimits, blockId)
  
  deallocate(xCoord)
  deallocate(yCoord)
  deallocate(zCoord)
  
  return

end subroutine Simulation_initBlock
