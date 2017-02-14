!!****if* source/Simulation/SimulationMain/magnetoHD/FieldLoop/Simulation_initBlock
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
  real :: x1,x2,x3,cos_ang,sin_ang,lambda
  real :: xx,yy,zz
#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_IHI_GC+1,GRID_JHI_GC+1,GRID_KHI_GC+1) :: Az,Ax,Ay
#else
  real, allocatable, dimension(:,:,:) :: Az,Ax,Ay
#endif

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

#ifndef FIXEDBLOCKSIZE
  if (NDIM == 2) then
     allocate(Ax(sizeX+1,sizeY+1,1),stat=istat)
     allocate(Ay(sizeX+1,sizeY+1,1),stat=istat)
     allocate(Az(sizeX+1,sizeY+1,1),stat=istat)
  elseif (NDIM == 3) then
     allocate(Ax(sizeX+1,sizeY+1,sizeZ+1),stat=istat)
     allocate(Ay(sizeX+1,sizeY+1,sizeZ+1),stat=istat)
     allocate(Az(sizeX+1,sizeY+1,sizeZ+1),stat=istat)
  endif
#endif

  if (NDIM == 3) then
     call Grid_getCellCoords(KAXIS,blockId,CENTER,    sim_gCell,zCoord, sizeZ)
     call Grid_getCellCoords(KAXIS,blockId,LEFT_EDGE, sim_gCell,zCoordL,sizeZ)
     call Grid_getCellCoords(KAXIS,blockId,RIGHT_EDGE,sim_gCell,zCoordR,sizeZ)
  endif
  if (NDIM >= 2) then
     call Grid_getCellCoords(JAXIS,blockId,CENTER,    sim_gCell,yCoord, sizeY)
     call Grid_getCellCoords(JAXIS,blockId,LEFT_EDGE, sim_gCell,yCoordL,sizeY)
     call Grid_getCellCoords(JAXIS,blockId,RIGHT_EDGE,sim_gCell,yCoordR,sizeY)
  endif

  call Grid_getCellCoords(IAXIS,blockId,CENTER,    sim_gCell,xCoord, sizeX)
  call Grid_getCellCoords(IAXIS,blockId,LEFT_EDGE, sim_gCell,xCoordL,sizeX)
  call Grid_getCellCoords(IAXIS,blockId,RIGHT_EDGE,sim_gCell,xCoordR,sizeX)

  !------------------------------------------------------------------------------
  ! Construct Az at each cell corner
  ! Bx = dAz/dy - dAy/dz
  ! By = dAx/dz - dAz/dx
  ! Bz = dAy/dx - dAx/dy
  Az = 0.
  Ax = 0.
  Ay = 0.

  rot = atan(sim_rx/sim_ry)
  cos_ang = cos(rot)
  sin_ang = sin(rot)

  if (cos_ang >= sin_ang) then
     lambda = (sim_xMax-sim_xMin)*cos_ang
  else
     lambda = (sim_zMax-sim_zMin)*sin_ang
  endif


  call Grid_getDeltas(blockID,del)
  dx = del(1)
  dy = del(2)
  dz = del(3)


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
           ! define radius of the field loop
           radius = sqrt((x1-sim_xCtr)**2 + (x2-sim_yCtr)**2)

           if (radius <= sim_fieldLoopRadius) then
              Ax(i,j,k) = 0.
              Ay(i,j,k) = 0.
              Az(i,j,k) = sim_Az_initial*(sim_fieldLoopRadius - radius)
           else
              Ax(i,j,k) = 0.
              Ay(i,j,k) = 0.
              Az(i,j,k) = 0.
           endif
        enddo
     enddo
  elseif (NDIM == 3) then
     do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)+1
        do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)+1
           do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)+1

              !! Rotated coordinate system
              !! / x1 \    / cos  0  sin \
              !! | x2 |  = |  0   1   0  |
              !! \ x3 /    \ -sin 0  cos /

#if NFACE_VARS > 0
              ! x Coord at cell corner
              if (i <=blkLimitsGC(HIGH,IAXIS)) then
                 xx = xCoordL(i)
              else
                 xx = xCoordR(i-1)
              endif

              ! y Coord at cell corner
              if (j <=blkLimitsGC(HIGH,JAXIS)) then
                 yy = yCoordL(j)
              else
                 yy = yCoordR(j-1)
              endif

              ! z Coord at cell corner
              if (k <=blkLimitsGC(HIGH,KAXIS)) then
                 zz = zCoordL(k)
              else
                 zz = zCoordR(k-1)
              endif

#else
              ! x Coord at cell center
              if (i <=blkLimitsGC(HIGH,IAXIS)) then
                 xx = xCoord(i)
              else
                 xx = xCoord(i-1) + dx
              endif

              ! y Coord at cell center
              if (j <=blkLimitsGC(HIGH,JAXIS)) then
                 yy = yCoord(j)
              else
                 yy = yCoord(j-1) + dy
              endif

              ! z Coord at cell center
              if (k <=blkLimitsGC(HIGH,KAXIS)) then
                 zz = zCoord(k)
              else
                 zz = zCoord(k-1) + dz
              endif
#endif
              ! For Ax and Ay
              x1 = (cos_ang*xCoord(i) + sin_ang*zz) ! with rotation
              !x1 = xx !without any rotation
              x2 = yy


              do while (x1 > 0.5*lambda)
                 x1 = x1 - lambda
              enddo
              do while (x1 < -0.5*lambda)
                 x1 = x1 + lambda
              enddo

              radius = sqrt(x1**2 + x2**2)

              if (radius < sim_fieldLoopRadius) then
                 Ax(i,j,k) = sim_Az_initial*(sim_fieldLoopRadius - radius)*(-sin_ang)
              endif
              Ay(i,j,k) = 0.

              ! For Az
              x1 = (cos_ang*xx + sin_ang*zCoord(k)) ! with rotation
              !x1 = xx !without any rotation
              x2 = yy


              do while (x1 > 0.5*lambda)
                 x1 = x1 - lambda
              enddo
              do while (x1 < -0.5*lambda)
                 x1 = x1 + lambda
              enddo

              radius = sqrt(x1**2 + x2**2)

              if (radius < sim_fieldLoopRadius) then
                 Az(i,j,k) = sim_Az_initial*(sim_fieldLoopRadius - radius)*cos_ang
              endif




           enddo
        enddo
     enddo
  endif



  call Grid_getBlkPtr(blockID,solnData,CENTER)

#if NFACE_VARS > 0
  if (sim_killdivb) then
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
     if (NDIM == 3) call Grid_getBlkPtr(blockID,facezData,FACEZ)
  endif
#endif


  do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)

           solnData(DENS_VAR,i,j,k)=  1.
           if (NDIM == 2) then
              solnData(VELX_VAR,i,j,k)= sim_U_initial*cos_ang
              solnData(VELY_VAR,i,j,k)= sim_U_initial*sin_ang
              solnData(VELZ_VAR,i,j,k)= sim_velz_initial
           elseif (NDIM == 3) then
              ! original version
              solnData(VELX_VAR,i,j,k)= sim_U_initial    ! 1.0 is a default.
              solnData(VELY_VAR,i,j,k)= sim_U_initial    ! 1.0 is a default.
              solnData(VELZ_VAR,i,j,k)= sim_velz_initial ! 2.0 is a default.
!!$
!!$              ! Dongwook's modification
!!$              solnData(VELX_VAR,i,j,k)= sim_U_initial*cos_ang    ! 1.0 is a default.
!!$              solnData(VELY_VAR,i,j,k)= sim_U_initial*sin_ang    ! 1.0 is a default.
!!$              solnData(VELZ_VAR,i,j,k)= sim_velz_initial ! 2.0 is a default.

           endif

           solnData(PRES_VAR,i,j,k)=  1.
           solnData(TEMP_VAR,i,j,k)=  1.

#if NFACE_VARS == 0
           solnData(MAGX_VAR,i,j,k)=  .5*(Az(i,j+1,k)-Az(i,j,k) + Az(i+1,j+1,k)-Az(i+1,j,k))/dy
           solnData(MAGY_VAR,i,j,k)= -.5*(Az(i+1,j,k)-Az(i,j,k) + Az(i+1,j+1,k)-Az(i,j+1,k))/dx
           solnData(MAGZ_VAR,i,j,k)= 0.
           solnData(MAGP_VAR,i,j,k) = .5*dot_product(solnData(MAGX_VAR:MAGZ_VAR,i,j,k),&
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
                 facexData(MAG_FACE_VAR,i,j,k)= (Az(i,j+1,k)-Az(i,j,k))/dy
                 faceyData(MAG_FACE_VAR,i,j,k)=-(Az(i+1,j,k)-Az(i,j,k))/dx
              elseif (NDIM == 3) then
                 facexData(MAG_FACE_VAR,i,j,k)= -(Ay(i,j,k+1)-Ay(i,j,k))/dz + (Az(i,j+1,k)-Az(i,j,k))/dy
                 faceyData(MAG_FACE_VAR,i,j,k)=  (Ax(i,j,k+1)-Ax(i,j,k))/dz - (Az(i+1,j,k)-Az(i,j,k))/dx
                 facezData(MAG_FACE_VAR,i,j,k)= -(Ax(i,j+1,k)-Ax(i,j,k))/dy + (Ay(i+1,j,k)-Ay(i,j,k))/dx
              endif
           endif
#else
           !! In this case we initialized Az using the cell-centered coordinates.
           if (NDIM == 2) then
              solnData(MAGX_VAR,i,j,k)= 0.5*(Az(i,j+1,k)-Az(i,j-1,k))/dy
              solnData(MAGY_VAR,i,j,k)=-0.5*(Az(i+1,j,k)-Az(i-1,j,k))/dx
           elseif (NDIM == 3) then
              solnData(MAGX_VAR,i,j,k)= -0.5*((Ay(i,j,k+1)-Ay(i,j,k-1))/dz + (Az(i,j+1,k)-Az(i,j-1,k))/dy)
              solnData(MAGY_VAR,i,j,k)=  0.5*((Ax(i,j,k+1)-Ax(i,j,k-1))/dz - (Az(i+1,j,k)-Az(i-1,j,k))/dx)
              solnData(MAGZ_VAR,i,j,k)= -0.5*((Ax(i,j+1,k)-Ax(i,j-1,k))/dy + (Ay(i+1,j,k)-Ay(i-1,j,k))/dx)
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

#if NDIM == 1
           solnData(DIVB_VAR,i,j,k) = 0.
#elif NDIM >= 2
           solnData(DIVB_VAR,i,j,k)= &
                     (facexData(MAG_FACE_VAR,i+1,j,  k  ) - facexData(MAG_FACE_VAR,i,j,k))/dx &
                   + (faceyData(MAG_FACE_VAR,i,  j+1,k  ) - faceyData(MAG_FACE_VAR,i,j,k))/dy
#if NDIM == 3
           solnData(DIVB_VAR,i,j,k)= solnData(DIVB_VAR,i,j,k) &
                   + (facezData(MAG_FACE_VAR,i,  j,  k+1) - facezData(MAG_FACE_VAR,i,j,k))/dz
#endif
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


#ifdef CURJ_VAR
  k=1
  do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
     do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
        solnData(CURJ_VAR,i,j,k) =&
               (faceyData(MAG_FACE_VAR,i+1,j,k)-faceyData(MAG_FACE_VAR,i,j,k))/dx &
              -(facexData(MAG_FACE_VAR,i,j+1,k)-facexData(MAG_FACE_VAR,i,j,k))/dy
     enddo
  enddo
#endif

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



