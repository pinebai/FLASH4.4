!!****if* source/Simulation/SimulationMain/unitTest/Gravity/BHTree-jeans/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!! 
!! SYNOPSIS
!!
!!  xall Simulation_initBlock(integer :: blockId)
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.  This version sets up the spherical expanding shell.
!!
!!
!! ARGUMENTS
!!
!!  blockId -        The number of the block to initialize
!!
!! PARAMETERS
!!
!!
!!
!!***

subroutine Simulation_initBlock (blockId)

  use Simulation_data
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_releaseBlkPtr,&
    Grid_getCellCoords, Grid_putPointData
  
  implicit none

#include "constants.h"
#include "Flash.h"
  
  integer,intent(IN) ::  blockId

  
  integer  :: i, j, k, n
  integer  :: ii, jj, kk, jhi, jlo
  real     :: vx, vy, vz, p, rho, e, ek, gpot
  real     :: xx, yy, zz, kx, k2

  real,allocatable,dimension(:) :: xCoord,yCoord,zCoord
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
  integer :: sizeX,sizeY,sizeZ
  integer,dimension(MDIM) :: axis
  real, dimension(:,:,:,:),pointer :: solnData
  integer,dimension(MDIM) :: startingPos
  real          :: del(MDIM)

  logical :: gcell = .true.

     
!! ---------------------------------------------------------------------------

  ! get the coordinate information for the current block from the database
  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
  sizeX = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
  allocate(xCoord(sizeX))
  sizeY = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
  allocate(yCoord(sizeY))
  sizeZ = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1
  allocate(zCoord(sizeZ))

  if (NDIM == 3) call Grid_getCellCoords&
                      (KAXIS, blockId, CENTER, gcell, zCoord, sizeZ)
  if (NDIM >= 2) call Grid_getCellCoords&
                      (JAXIS, blockId, CENTER,gcell, yCoord, sizeY)
  call Grid_getCellCoords(IAXIS, blockId, CENTER, gcell, xCoord, sizeX)
  call Grid_getDeltas(blockId,del)
  !
  !     For each cell
  !  
#ifdef FL_NON_PERMANENT_GUARDCELLS
  call Grid_getBlkPtr(blockId,solnData)
#endif


  do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
    zz = zCoord(k)
    do j = blkLimitsGC(LOW, JAXIS), blkLimitsGC(HIGH, JAXIS)
       yy = yCoord(j)
       do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH, IAXIS)
         xx = xCoord(i)
             
         kx   = 2*PI*(sim_hx*xx/sim_Lx + sim_hy*yy/sim_Ly + sim_hz*zz/sim_Lz)
         k2   = 4*PI*PI*(sim_hx*sim_hx/sim_Lx/sim_Lx &
         &    + sim_hy*sim_hy/sim_Ly/sim_Ly + sim_hz*sim_hz/sim_Lz/sim_Lz)
         rho  = sim_rho0*(1+sim_delta*cos(kx))
         p    = sim_p0  *(1+sim_delta*cos(kx))
         vx   = 0.0
         vy   = 0.0
         vz   = 0.0
         gpot = -4*PI*sim_newton*sim_rho0*sim_delta/k2 * cos(kx)


         ! assume gamma-law equation of state
         ek  = 0.5*(vx*vx + vy*vy + vz*vz)
         e   = p/(sim_gamma_1-1.)
         e   = e/rho + ek
         e   = max (e, smallP)
           
         axis(IAXIS)=i
         axis(JAXIS)=j
         axis(KAXIS)=k


#ifdef FL_NON_PERMANENT_GUARDCELLS
         solnData(DENS_VAR,i,j,k)=rho
         solnData(PRES_VAR,i,j,k)=p
         solnData(ENER_VAR,i,j,k)=e
         solnData(GAME_VAR,i,j,k)=sim_gamma_1
         solnData(GAMC_VAR,i,j,k)=sim_gamma_1
         solnData(VELX_VAR,i,j,k)=vx
         solnData(VELY_VAR,i,j,k)=vy
         solnData(VELZ_VAR,i,j,k)=vz
         solnData(PANL_VAR,i,j,k)=gpot
#else

         call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rho)
         call Grid_putPointData(blockId, CENTER, PRES_VAR, EXTERIOR, axis, p)
         call Grid_putPointData(blockId, CENTER, ENER_VAR, EXTERIOR, axis, e)    
         call Grid_putPointData(blockId, CENTER, GAME_VAR, EXTERIOR, axis, sim_gamma_1)
         call Grid_putPointData(blockId, CENTER, GAMC_VAR, EXTERIOR, axis, sim_gamma_1)
         call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, axis, vx)
         call Grid_putPointData(blockId, CENTER, VELY_VAR, EXTERIOR, axis, vy)
         call Grid_putPointData(blockId, CENTER, VELZ_VAR, EXTERIOR, axis, vz)
         call Grid_putPointData(blockId, CENTER, PANL_VAR, EXTERIOR, axis, gpot)
#endif
        enddo
     enddo
  enddo
#ifdef FL_NON_PERMANENT_GUARDCELLS
  call Grid_releaseBlkPtr(blockID, solnData)
#endif
  deallocate(xCoord)
  deallocate(yCoord)
  deallocate(zCoord)
  return
end subroutine Simulation_initBlock






