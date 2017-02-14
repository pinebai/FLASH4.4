!!****if* source/Simulation/SimulationMain/SinkRotatingCloudCore/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!! SYNOPSIS
!!
!!  call Simulation_initBlock(integer,intent(IN)  :: blockid)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   blockid : ID of block in current processor
!!
!! AUTOGENROBODOC
!!
!!
!!***



subroutine Simulation_initBlock (blockId)

  use Simulation_data, ONLY : sim_xMin, sim_yMin, sim_zMin, &
       sim_xMax, sim_yMax, sim_zMax
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
       Grid_getCellCoords, Grid_putPointData
  use Driver_interface, ONLY : Driver_abortFlash
  use Eos_interface, ONLY : Eos_wrapped
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get  

  implicit none
  
#include "constants.h"
#include "Flash.h"
#include "Eos.h"
  
  integer,intent(IN) ::  blockId
  
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer :: sizeX, sizeY, sizeZ
  integer, dimension(MDIM) :: axis
  real, dimension(:), allocatable :: x, y, z
  logical :: gcell=.true.
  integer :: i,j,k,istat
   
  real :: radius, xcenter, ycenter, zcenter, radiusSphere
  real :: rho, temp, vx, vy, vz, m
  real :: xdist, ydist, zdist, distxy
  real :: phi, pres, cs2
  real :: bb_radius, bb_dens, bb_cs, bb_omega

  call RuntimeParameters_get("bb_radius", bb_radius)
  call RuntimeParameters_get("bb_dens", bb_dens)  
  call RuntimeParameters_get("bb_cs", bb_cs)
  call RuntimeParameters_get("bb_omega", bb_omega)

  call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)
  sizeX = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS) + 1
  sizeY = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS) + 1
  sizeZ = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS) + 1
  allocate(x(sizeX), stat=istat)
  allocate(y(sizeY), stat=istat)
  allocate(z(sizeZ), stat=istat)
  call Grid_getCellCoords(KAXIS, blockId, CENTER, gcell, z, sizeZ)
  call Grid_getCellCoords(JAXIS, blockId, CENTER, gcell, y, sizeY)
  call Grid_getCellCoords(IAXIS, blockId, CENTER, gcell, x, sizeX)

  xcenter = 0.0
  ycenter = 0.0
  zcenter = 0.0
  
  do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
     do j = blkLimitsGC(LOW, JAXIS), blkLimitsGC(HIGH, JAXIS)
        do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)

           xdist = x(i) - xcenter
           ydist = y(j) - ycenter
           zdist = z(k) - zcenter
    
           radius = sqrt( xdist**2 + ydist**2 + zdist**2)
           distxy =  sqrt( xdist**2 + ydist**2 )
                      
           cs2 = bb_cs
           cs2 = cs2*cs2
           
           if(xdist .ne. 0.0) then
              phi = ATAN(ydist/xdist)
           else
              phi = 3.1415926535/2.0
           end if
           
           if(radius .le. bb_radius) then
              rho = bb_dens
              ! m=2 10% perturbation:
              rho = rho * (1.0 + 0.1 * cos(2.0 * phi))
              vx = -distxy * bb_omega * ABS(SIN(phi)) * &
                   SIGN(1.0, ydist)
              vy = distxy * bb_omega * COS(phi) * &
                   SIGN(1.0, xdist)
              pres = cs2*rho
           else
              rho = 1e-2*bb_dens
              vx = 0.0
              vy = 0.0
              pres = cs2*bb_dens
           end if

           vz = 0.0
    
           axis(IAXIS) = i
           axis(JAXIS) = j
           axis(KAXIS) = k
           
           call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rho)
           call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, axis, vx)
           call Grid_putPointData(blockId, CENTER, VELY_VAR, EXTERIOR, axis, vy)
           call Grid_putPointData(blockId, CENTER, VELZ_VAR, EXTERIOR, axis, vz)
           call Grid_putPointData(blockId, CENTER, PRES_VAR, EXTERIOR, axis, pres)

        end do
     end do
  end do
  
  call Eos_wrapped(MODE_DENS_PRES, blklimitsGC, blockId)

  deallocate(x)
  deallocate(y)
  deallocate(z)

  return

end subroutine Simulation_initBlock
