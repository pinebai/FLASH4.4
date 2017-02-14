!!****if* source/Simulation/SimulationMain/Jeans/Simulation_initBlock
!!
!! NAME
!! 
!! Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer :: blockId,
!!                       
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.  This version sets up the Jeans instability 
!!  problem.
!!
!!  Reference:   Chandrasekhar, S., 1961, Hydrodynamic and Hydromagnetic
!!                  Stability (Oxford:  OUP), p. 588
!!
!! ARGUMENTS
!!
!!  blockId -        The number of the block to initialize
!!
!!***
subroutine Simulation_initBlock(blockID)

  use Simulation_data, ONLY: sim_kZ, sim_kY, sim_kX, sim_A, sim_velA, &
       &  sim_rho0, sim_p0, sim_smallE, sim_gamma
  use Grid_interface, ONLY: Grid_getBlkIndexLimits, Grid_putRowData, &
       Grid_getCellCoords

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent(IN)  ::  blockID

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer, dimension(3) :: startingPos
  integer :: i, j, k, n
  integer :: sizeX, sizeY, sizeZ

  logical :: gcell = .true.

  real    ::  xcos, ycos, zcos, xsin, ysin, zsin
  real    ::  xarg, yarg, zarg
  real    ::  pert, vpert
  
  real, dimension(:), allocatable :: x, y, z
  real, dimension(:), allocatable :: rho, p, t, e, ek, gam, vx, vy, vz, ei

  ! Get the coordinate information for the current block
      
  call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)
  sizeX = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS) + 1
  allocate(x(sizeX))
  sizeY = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS) + 1
  allocate(y(sizeY))
  sizeZ = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS) + 1
  allocate(z(sizeZ))
  if (NDIM == 3) then 
     call Grid_getCellCoords(KAXIS, blockId, CENTER, gcell, z, sizeZ)
  endif
  if (NDIM >= 2) then 
     call Grid_getCellCoords(JAXIS, blockId, CENTER, gcell, y, sizeY)
  endif
  call Grid_getCellCoords(IAXIS, blockId, CENTER, gcell, x, sizeX)

  ! Set initial conditions in each zone

  allocate(vx(sizeX), vy(sizeX), vz(sizeX), p(sizeX), rho(sizeX), e(sizeX), & 
       ek(sizeX), ei(sizeX), gam(sizeX))
  
  do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
     zarg = z(k)*sim_kZ*K3D
     zsin = sin(zarg)
     zcos = cos(zarg)
     if (NDIM < 3) then
        zsin = 0.
        zcos = 1.
     endif
     do j = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
        yarg = y(j)*sim_kY*K2D
        ysin = sin(yarg)
        ycos = cos(yarg)
        if (NDIM < 2) then
           ysin = 0.
           ycos = 1.
        endif
        do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)
           xarg = x(i)*sim_kX
           xsin = sin(xarg)
           xcos = cos(xarg)
           pert = sim_A * (zcos*(xcos*ycos+xsin*ysin) - & 
                zsin*(xcos*ysin+xsin*ycos))
           vpert = sim_velA * (zcos*(xcos*ysin+xsin*ycos) + & 
                zsin*(xcos*ycos+xsin*ysin))
           rho(i) = sim_rho0*(1.+pert)
           p  (i) = sim_p0*(1.+sim_gamma*pert)
           vx (i) = 0.
           vy (i) = 0.
           vz (i) = 0.

        enddo

        gam = sim_gamma
        ek = 0.5 * (vx*vx+vy*vy+vz*vz)
        ei  = p/(rho*(sim_gamma-1.))
        e  = ei + ek
        e  = max(e, sim_smallE)

        startingPos(1) = 1
        startingPos(2) = j
        startingPos(3) = k
        call Grid_putRowData(blockID, CENTER, DENS_VAR, EXTERIOR, &
             IAXIS, startingPos, rho, sizeX)
        call Grid_putRowData(blockID, CENTER, PRES_VAR, EXTERIOR, &
             IAXIS, startingPos, p, sizeX)
        call Grid_putRowData(blockID, CENTER, ENER_VAR, EXTERIOR, &
             IAXIS, startingPos, e, sizeX)
        call Grid_putRowData(blockID, CENTER, GAME_VAR, EXTERIOR, &
             IAXIS, startingPos, gam, sizeX)
        call Grid_putRowData(blockID, CENTER, GAMC_VAR, EXTERIOR, &
             IAXIS, startingPos, gam, sizeX)
        call Grid_putRowData(blockID, CENTER, VELX_VAR, EXTERIOR, &
             IAXIS, startingPos, vx, sizeX)
        call Grid_putRowData(blockID, CENTER, VELY_VAR, EXTERIOR, &
             IAXIS, startingPos, vy, sizeX)
        call Grid_putRowData(blockID, CENTER, VELZ_VAR, EXTERIOR, &
             IAXIS, startingPos, vz, sizeX)
        call Grid_putRowData(blockID, CENTER, EINT_VAR, EXTERIOR, &
             IAXIS, startingPos, ei, sizeX)
 
     enddo
  enddo

  deallocate(rho, p, vx, vy, vz, e, ei, ek, gam)
  deallocate(x)
  deallocate(y)
  deallocate(z)
  
  return
end subroutine Simulation_initBlock

