!!****if* source/Simulation/SimulationMain/magnetoHD/BeerCan/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(in) :: blockID) 
!!                       
!!
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.
!!
!!  Reference:  Gardiner & Stone JCP 205(2005),509-539
!!
!!  Parameters:  blockID      The number of the block to initialize
!!
!! 
!! ARGUMENTS
!!
!!  blockID -          the number of the block to update
!!
!! 
!!
!!***

subroutine Simulation_initBlock(blockID)

  use Simulation_data

  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
                             Grid_getCellCoords,     &
                             Grid_getDeltas, &
                             Grid_getBlkPtr, &
                             Grid_releaseBlkPtr
  implicit none

#include "constants.h"
#include "Flash.h"

  !! Arguments ------------------------
  integer, intent(in) :: blockID
  !! ----------------------------------

  integer :: i, j, k, n, istat, sizeX, sizeY, sizeZ
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real :: enerZone, ekinZone, eintZone, rot, radius, dx, dy, dz
  real, allocatable,dimension(:) :: xCoord,xCoordL,xCoordR,&
                                    yCoord,yCoordL,yCoordR,&
                                    zCoord,zCoordL,zCoordR
  real, dimension(MDIM) :: del
  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData
  real :: x1,x2,x3,cos_ang,sin_ang,lambda, Btor, f, g, c_light, cur
  real :: xx,yy,zz
#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_IHI_GC+1,GRID_JHI_GC+1,GRID_KHI_GC+1) :: Az,Ax,Ay
#else
  real, allocatable, dimension(:,:,:) :: Az,Ax,Ay
#endif

  ! dump some output to stdout listing the paramters
!!$   if (sim_meshMe == MASTER_PE) then
!!$1    format (1X, 1P, 4(A7, E13.7, :, 1X))
!!$2    format (1X, 1P, 2(A7, E13.7, 1X), A7, I13)
!!$  endif
  
  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)

  sizeX = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  sizeY = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
  sizeZ = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

  allocate(xCoord(sizeX), stat=istat)
  allocate(xCoordL(sizeX),stat=istat)
  allocate(xCoordR(sizeX),stat=istat)

  allocate(yCoord(sizeY), stat=istat)
  allocate(yCoordL(sizeY),stat=istat)
  allocate(yCoordR(sizeY),stat=istat)

  allocate(zCoord(sizeZ), stat=istat)
  allocate(zCoordL(sizeZ),stat=istat)
  allocate(zCoordR(sizeZ),stat=istat)

  xCoord  = 0.0
  xCoordL = 0.0
  xCoordR = 0.0

  yCoord  = 0.0
  yCoordL = 0.0
  yCoordR = 0.0

  zCoord  = 0.0
  zCoordL = 0.0
  zCoordR = 0.0


  if (NDIM == 3) then
     call Grid_getCellCoords(KAXIS,blockId,CENTER,    sim_gCell,zCoord, sizeZ)
  endif
  if (NDIM >= 2) then
     call Grid_getCellCoords(JAXIS,blockId,CENTER,    sim_gCell,yCoord, sizeY)
  endif

  call Grid_getCellCoords(IAXIS,blockId,CENTER,    sim_gCell,xCoord, sizeX)

  call Grid_getDeltas(blockID,del)
  dx = del(1)
  dy = del(2)
  dz = del(3)
  call Grid_getBlkPtr(blockID,solnData,CENTER)

  c_light = 2.99e10
  
  if (NDIM == 3) then
     do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
        do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
           do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)

              xx = xCoord(i)
              yy = yCoord(j)
              zz = zCoord(k)
           
              radius = sqrt(xx**2 + yy**2)
              cos_ang = xx/radius
              sin_ang = yy/radius
              cur = sim_b*c_light*sim_r0/2.
              ! Set value of f
              if (radius < sim_r0) then
                 f= 2.*cur*radius/sim_r0**2/c_light
              else if (radius <= sim_r1) then
                 f= 2.*cur/radius/c_light
              else if (radius <= sim_r1+sim_delta) then
                 f = 2.*cur/sim_r1/c_light * (1.-(radius-sim_r1)/sim_delta)
              else
                 f = 0.0
              endif 

              ! Set value of g
              if (abs(zz) < sim_z0/2.) then
                 g= 1.
              else if (abs(zz) <= sim_z0/2. +sim_delta) then
                 g= 1.- (abs(zz)-sim_z0/2.)/sim_delta
              else
                 g = 0.0
              endif 
              
              Btor = f*g
              !convert to flash units 
              Btor = Btor/sqrt(4.*PI)
              
              solnData(MAGX_VAR,i,j,k)=  -Btor*sin_ang
              solnData(MAGY_VAR,i,j,k)=   Btor*cos_ang
              solnData(MAGZ_VAR,i,j,k)=   0.0
        
              solnData(VELX_VAR,i,j,k)= 0.0
              solnData(VELY_VAR,i,j,k)= 0.0
              solnData(VELZ_VAR,i,j,k)= 0.0
              solnData(PRES_VAR,i,j,k)=  1.
              solnData(TEMP_VAR,i,j,k)=  1.
              solnData(DENS_VAR,i,j,k)=  1.

              ! Compute the gas energy and set the gamma-values needed for the EOS
              ekinZone = 0.5 * dot_product(solnData(VELX_VAR:VELZ_VAR,i,j,k),&
                                        solnData(VELX_VAR:VELZ_VAR,i,j,k))

              ! specific internal energy
              eintZone = solnData(PRES_VAR,i,j,k)/(sim_gamma-1.)/solnData(DENS_VAR,i,j,k)

              ! total specific gas energy
              enerZone = eintZone + ekinZone

              ! Take a limit value
              enerZone = max(enerZone, sim_smallP)

              solnData(ENER_VAR,i,j,k)=enerZone
              solnData(EINT_VAR,i,j,k)=eintZone
              solnData(GAMC_VAR,i,j,k)=sim_gamma
              solnData(GAME_VAR,i,j,k)=sim_gamma

           enddo
        enddo
     enddo
  endif


  ! Release pointer
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)

#if NFACE_VARS > 0
  if (sim_killdivb) then
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
     if (NDIM == 3) call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
  endif
#endif

  deallocate(xCoord)
  deallocate(xCoordL)
  deallocate(xCoordR)

  deallocate(yCoord)
  deallocate(yCoordL)
  deallocate(yCoordR)

  deallocate(zCoord)
  deallocate(zCoordL)
  deallocate(zCoordR)

#ifndef FIXEDBLOCKSIZE
  deallocate(Az)
  deallocate(Ax)
  deallocate(Ay)
#endif

end subroutine Simulation_initBlock



