!!****if* source/Simulation/SimulationMain/magnetoHD/NohCylindrical/Simulation_initBlock
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
!!  Reference: Velikovich et al. Phys. Plasmas 19 (2012), 012707
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
  real :: radius, dx, dy, dz
  real, allocatable,dimension(:) :: xCoord,xCoordL,xCoordR,&
                                    yCoord,yCoordL,yCoordR,&
                                    zCoord,zCoordL,zCoordR
  real :: enerZone, ekinZone, eintZone
  real, dimension(MDIM) :: del
  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData
  real :: x1,x2,x3
  real :: xx,yy,zz

#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_IHI_GC+1,GRID_JHI_GC+1,GRID_KHI_GC+1) :: Az,Ax,Ay
#else
  real, allocatable, dimension(:,:,:) :: Az,Ax,Ay
#endif

  
  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)

!* *********************************** *
!*  Get block limits to create grid    *
!* *********************************** *

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

!* *********************************** *
!*  Allocate vector potential arrays   *
!* *********************************** *

#ifndef FIXEDBLOCKSIZE
  if (NDIM == 2) then
     allocate(Ax(sizeX+1,sizeY+1,1),stat=istat)
     allocate(Ay(sizeX+1,sizeY+1,1),stat=istat)
     allocate(Az(sizeX+1,sizeY+1,1),stat=istat)
  endif
#endif

  Az = 0.
  Ax = 0.
  Ay = 0.

!* ******************************************* *
!*  get goordinates for cells (center/edges)   *
!* ******************************************* *

  if (NDIM >= 2) then
     call Grid_getCellCoords(JAXIS,blockId,CENTER,    sim_gCell,yCoord, sizeY)
     call Grid_getCellCoords(JAXIS,blockId,LEFT_EDGE, sim_gCell,yCoordL,sizeY)
     call Grid_getCellCoords(JAXIS,blockId,RIGHT_EDGE,sim_gCell,yCoordR,sizeY)
  endif

  call Grid_getCellCoords(IAXIS,blockId,CENTER,    sim_gCell,xCoord, sizeX)
  call Grid_getCellCoords(IAXIS,blockId,LEFT_EDGE, sim_gCell,xCoordL,sizeX)
  call Grid_getCellCoords(IAXIS,blockId,RIGHT_EDGE,sim_gCell,xCoordR,sizeX)

  call Grid_getDeltas(blockID,del)
  dx = del(1)
  dy = del(2)
  dz = del(3)

!* ************************** *
!* get coords at edges        *
!* ************************** *

  if (NDIM == 2) then
     k = 1
     do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)+1
        do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)+1

#if NFACE_VARS > 1
           ! x Coord at cell corner
           if (i <=blkLimitsGC(HIGH,IAXIS)) then
              x1 = xCoordL(i)
           else
              x1 = xCoordR(i-1)
           endif

           ! y Coord at cell corner
           if (j <=blkLimitsGC(HIGH,JAXIS)) then
              x2 = yCoordL(j)
           else
              x2 = yCoordR(j-1)
           endif
#else
           ! x Coord at cell center
           if (i <=blkLimitsGC(HIGH,IAXIS)) then
              x1 = xCoord(i)
           else
              x1 = xCoord(i-1) + dx
           endif

           ! y Coord at cell center
           if (j <=blkLimitsGC(HIGH,JAXIS)) then
              x2 = yCoord(j)
           else
              x2 = yCoord(j-1) + dy
           endif
#endif


           ! define radius with respect to center
           radius = x1
           ! define Az  
           Ay(i,j,k) =  -0.5*radius**2*6.35584*1.e5/sqrt(4.*PI*sim_UnitDensity*sim_UnitVelocity**2)

        enddo
     enddo
  endif

  call Grid_getBlkPtr(blockID,solnData,CENTER)

#if NFACE_VARS > 0
  if (sim_killdivb) then
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
  endif
#endif

!* ************************************ *
! INIT, set primitives  
!
!* ************************************ *

  do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)

           ! x Coord at cell center
           if (i <=blkLimitsGC(HIGH,IAXIS)) then
              x1 = xCoord(i)
           else
              x1 = xCoord(i-1) + dx
           endif

           ! y Coord at cell center
           if (j <=blkLimitsGC(HIGH,JAXIS)) then
              x2 = yCoord(j)
           else
              x2 = yCoord(j-1) + dy
           endif

           !radius 
           radius = x1

           solnData(DENS_VAR,i,j,k)=  3.1831*radius**2
           if (NDIM == 2) then
              solnData(VELX_VAR,i,j,k)= -3.24101
              solnData(VELY_VAR,i,j,k)= 0.0
              solnData(VELZ_VAR,i,j,k)= 0.0
           endif     
         ! plasma beta is 8 PI 10^-6
           solnData(PRES_VAR,i,j,k)=  8.*PI*1.e-6*&
             (x1*6.35584*1.e5/sqrt(4.*PI*sim_UnitDensity*sim_UnitVelocity**2))**2

           solnData(TEMP_VAR,i,j,k)=  solnData(PRES_VAR,i,j,k)/solnData(DENS_VAR,i,j,k)


#if NFACE_VARS == 0
           solnData(MAGX_VAR,i,j,k)=  0.0
           solnData(MAGY_VAR,i,j,k)=  0.0
           solnData(MAGZ_VAR,i,j,k)=  x1*6.35584*1.e5/sqrt(4.*PI*sim_UnitDensity*sim_UnitVelocity**2)
           solnData(MAGP_VAR,i,j,k) = 0.5*dot_product(solnData(MAGX_VAR:MAGZ_VAR,i,j,k),&
                                                      solnData(MAGX_VAR:MAGZ_VAR,i,j,k))
           solnData(DIVB_VAR,i,j,k) = 0.

#endif
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

#ifdef FLASH_3T
           !Split energy, set erad to a small value
           solnData(EELE_VAR,i,j,k) = eintZone / 2.0 
           solnData(EION_VAR,i,j,k) = eintZone / 2.0 
           solnData(ERAD_VAR,i,j,k) = eintZone * 1.0e-12
#endif

#ifdef VECZ_VAR
           ! vector potential Az
           if (NFACE_VARS > 1) then
              solnData(VECZ_VAR,i,j,k) = .25*(Az(i,j,k)+Az(i+1,j,k)+Az(i,j+1,k)+Az(i+1,j+1,k))
           else
              solnData(VECZ_VAR,i,j,k) = Az(i,j,k)
           endif
#endif

#if NFACE_VARS > 0
           !! In this case we initialized Az using the cell-cornered coordinates.

           if (sim_killdivb) then
              if (NDIM == 2) then
                 facexData(MAG_FACE_VAR,i,j,k)=0.0! (Az(i,j+1,k)-Az(i,j,k))/dy
                 faceyData(MAG_FACE_VAR,i,j,k)=0.0!-(Az(i+1,j,k)-Az(i,j,k))/dx
                 solnData(MAGX_VAR,i,j,k)=  0.0
                 solnData(MAGY_VAR,i,j,k)=  0.0
                 solnData(MAGZ_VAR,i,j,k)=  x1*6.35584*1.e5/sqrt(4.*PI*sim_UnitDensity*sim_UnitVelocity**2)
              endif
           endif
#else
           !! In this case we initialized Az using the cell-centered coordinates.

           if (NDIM == 2) then
              solnData(MAGX_VAR,i,j,k)= 0.5*(Az(i,j+1,k)-Az(i,j-1,k))/dy
              solnData(MAGY_VAR,i,j,k)=-0.5*(Az(i+1,j,k)-Az(i-1,j,k))/dx
           endif
#endif

        enddo
     enddo
  enddo



  do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

#if NFACE_VARS > 0
           solnData(MAGX_VAR,i,j,k) = 0.5*(facexData(MAG_FACE_VAR,i,j,k)+facexData(MAG_FACE_VAR,i+1,j,k))
           solnData(MAGY_VAR,i,j,k) = 0.5*(faceyData(MAG_FACE_VAR,i,j,k)+faceyData(MAG_FACE_VAR,i,j+1,k))
           if (NDIM == 3) then
              solnData(MAGZ_VAR,i,j,k) = 0.5*(facezData(MAG_FACE_VAR,i,j,k)+facezData(MAG_FACE_VAR,i,j,k+1))
           endif

! initialize divB to zero (as B=B_phi only there is no stag. field to get something different...)           
#if NDIM == 1
           solnData(DIVB_VAR,i,j,k) = 0.0
#elif NDIM >= 2
           solnData(DIVB_VAR,i,j,k)=  0.0
#endif

#else !NFACE_VARS == 0
           solnData(DIVB_VAR,i,j,k) = 0.
#endif !NFACE_VARS

           ! Update the magnetic pressure
           solnData(MAGP_VAR,i,j,k) = .5*dot_product(solnData(MAGX_VAR:MAGZ_VAR,i,j,k),&
                                                     solnData(MAGX_VAR:MAGZ_VAR,i,j,k))

        enddo
     enddo
  enddo



  ! Release pointer
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)

#if NFACE_VARS > 0
  if (sim_killdivb) then
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
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



