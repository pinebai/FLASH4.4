!!****if* source/Simulation/SimulationMain/magnetoHD/Torus/Simulation_initBlock
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

  integer :: i, j, k, n, sizeX, sizeY, sizeZ
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real :: enerZone, ekinZone, eintZone, rot, radius, dx, dy, dz
  real, allocatable,dimension(:) :: xCoord,xCoordL,xCoordR,&
                                    yCoord,yCoordL,yCoordR,&
                                    zCoord,zCoordL,zCoordR
  real, dimension(MDIM) :: del
  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData
  real :: x1,x2,x3,cos_ang,sin_ang,lambda
  real :: xx,yy,zz

  real :: R, z, A_phi, gmm1
  real :: a,b,c, p0, rho0, kappa, l_k, l_k2, r_maxpr
  real :: Bo, epsilon, beta, Vo, R_in, rho_max, a_max, Tt_max
  real :: prt, pra, sch, phi, Temp_contr, rhoa_one, integral_g_spline
  real :: rhot, rhoa, ETA
  real :: R_sphere
  real :: ag,bg,cg
  real :: rhocut



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
  
  
  gmm1 = sim_gamma- 1.0

  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)

  sizeX = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  sizeY = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
  sizeZ = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

  allocate(xCoord(sizeX))
  allocate(xCoordL(sizeX))
  allocate(xCoordR(sizeX))

  allocate(yCoord(sizeY))
  allocate(yCoordL(sizeY))
  allocate(yCoordR(sizeY))

  allocate(zCoord(sizeZ))
  allocate(zCoordL(sizeZ))
  allocate(zCoordR(sizeZ))

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
     allocate(Ax(sizeX+1,sizeY+1,1))
     allocate(Ay(sizeX+1,sizeY+1,1))
     allocate(Az(sizeX+1,sizeY+1,1))
  elseif (NDIM == 3) then
     allocate(Ax(sizeX+1,sizeY+1,sizeZ+1))
     allocate(Ay(sizeX+1,sizeY+1,sizeZ+1))
     allocate(Az(sizeX+1,sizeY+1,sizeZ+1))
  endif
#endif

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



           radius = abs(x1)
           z = x2
           R = sqrt(radius**2 + z**2)
           phi = 0.0
           R_sphere = sim_rSphere
           r_maxpr = sim_rMax
           R_in = sim_rMin
           l_k = sqrt(r_maxpr)*r_maxpr/(r_maxpr-sim_r0)
           l_k2 = l_k*l_k
           Vo = l_k
           rhocut = sim_denCut
        
   ! *************************************
   ! Solve for kappa and c assuming pra=0
   !  
   ! ************************************* */  
 
           c =  -1.0/(R_in-sim_r0) + 0.5*l_k2/(R_in*R_in);
           kappa = (gmm1/sim_gamma)*(c + 0.5*l_k2*(-1./r_maxpr/r_maxpr) & 
                            + 1./(r_maxpr-sim_r0))/(sim_denMax**gmm1)
        
   ! ************************************
   ! Torus density + pressure
   ! ************************************* */  

           a          = c + 1.0/(R-sim_r0) - 0.5*l_k2/(radius**2)
           a          = max(a, 0.0)
           rhot       = (0.4*a/kappa)**1.5
           prt        = kappa*rhot**sim_gamma
           ETA        = sim_dCon
           Temp_contr = sim_tCon
           Tt_max     = kappa*sim_denMax**gmm1
        
   ! ************************************
   ! atmospheric density + pressure 
   ! (force equilibrium Fg and \nabla P)
   ! ************************************ */

           rhoa = ETA*sim_denMax*exp( ( 1./(R-sim_r0) - 1./(r_maxpr-sim_r0) )/Temp_contr/Tt_max ) 
           pra  = rhoa*Tt_max*Temp_contr
           rhoa_one = ETA*sim_denMax*exp( ( 1./(R_sphere-sim_r0) - 1./(r_maxpr-sim_r0) )/Temp_contr/Tt_max )

           rhoa = ETA*sim_denMax*exp( ( 1./(R-sim_r0) - 1./(R_in-sim_r0) )/Temp_contr/Tt_max)
           pra  = rhoa*Temp_contr*Tt_max
           Bo = sqrt(2.*kappa*(sim_denMax**sim_gamma)/sim_beta)/sim_denMax
           b = Vo/radius
           
           if (rhot > rhocut .AND. radius > 2.0) then           
             A_phi = Bo*(rhot - rhocut)             
           else
             A_phi = 0.0
           endif
           
           Ax(i,j,k) = 0.0
           Ay(i,j,k) = 0.0
           Az(i,j,k) = A_phi



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

!! INIT LOOP DEFINE STUFF HERE!



  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

           radius = abs(xCoord(i))
           z = yCoord(j)
           R = sqrt(radius**2 + z**2)
           phi = 0.0
           R_sphere = sim_rSphere
           r_maxpr = sim_rMax
           R_in = sim_rMin
           l_k = sqrt(r_maxpr)*r_maxpr/(r_maxpr-sim_r0)
           l_k2 = l_k*l_k
           Vo = l_k
           rhocut = sim_denCut
        
   ! *************************************
   ! Solve for kappa and c assuming pra=0
   !  
   ! ************************************* */  
 
           c =  -1.0/(R_in-sim_r0) + 0.5*l_k2/(R_in*R_in);
           kappa = (gmm1/sim_gamma)*(c + 0.5*l_k2*(-1./r_maxpr/r_maxpr) & 
                            + 1./(r_maxpr-sim_r0))/(sim_denMax**gmm1)
        
   ! ************************************
   ! Torus density + pressure
   ! ************************************* */  

           a          = c + 1.0/(R-sim_r0) - 0.5*l_k2/(radius**2)
           a          = max(a, 0.0)
           rhot       = (0.4*a/kappa)**1.5
           prt        = kappa*rhot**sim_gamma
           ETA        = sim_dCon
           Temp_contr = sim_tCon
           Tt_max     = kappa*sim_denMax**gmm1
        
   ! ************************************
   ! atmospheric density + pressure 
   ! (force equilibrium Fg and \nabla P)
   ! ************************************ */

           rhoa = ETA*sim_denMax*exp( ( 1./(R-sim_r0) - 1./(r_maxpr-sim_r0) )/Temp_contr/Tt_max ) 
           pra  = rhoa*Tt_max*Temp_contr
           rhoa_one = ETA*sim_denMax*exp( ( 1./(R_sphere-sim_r0) - 1./(r_maxpr-sim_r0) )/Temp_contr/Tt_max )

           rhoa = ETA*sim_denMax*exp( ( 1./(R-sim_r0) - 1./(R_in-sim_r0) )/Temp_contr/Tt_max)
           pra  = rhoa*Temp_contr*Tt_max
           Bo = sqrt(2.*kappa*(sim_denMax**sim_gamma)/sim_beta)/sim_denMax
           b = Vo/radius
           
           if (rhot > rhocut .AND. radius > 2.0) then           
             A_phi = Bo*(rhot - rhocut)
           else
             A_phi = 0.0
           endif
           
                      
           if (prt > pra .AND. radius > 2.0) then
              solnData(DENS_VAR,i,j,k)=  rhot
              solnData(VELX_VAR,i,j,k)=  0.0
              solnData(VELY_VAR,i,j,k)=  0.0
              solnData(VELZ_VAR,i,j,k)=  b
              solnData(PRES_VAR,i,j,k)=  prt
              solnData(TEMP_VAR,i,j,k)=  prt/rhot
           else
              solnData(DENS_VAR,i,j,k)=  rhoa
              solnData(VELX_VAR,i,j,k)=  0.0
              solnData(VELY_VAR,i,j,k)=  0.0
              solnData(VELZ_VAR,i,j,k)=  0.0
              solnData(PRES_VAR,i,j,k)=  pra
              solnData(TEMP_VAR,i,j,k)=  pra/rhoa
           endif 

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


  do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)

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
!                 facexData(MAG_FACE_VAR,i,j,k)= 1./(xCoordL(i)) !-(Az(i,j+1,k)-Az(i,j,k))/dy
!                 faceyData(MAG_FACE_VAR,i,j,k)= 1.0!(xCoordR(i)*Az(i+1,j,k)-xCoordL(i)*Az(i,j,k))/dx/xCoord(i)
                 facexData(MAG_FACE_VAR,i,j,k)= -(Az(i,j+1,k)-Az(i,j,k))/dy
                 faceyData(MAG_FACE_VAR,i,j,k)= (xCoordR(i)*Az(i+1,j,k)-xCoordL(i)*Az(i,j,k))/dx/xCoord(i)
!                 facexData(MAG_FACE_VAR,i,j,k)= sin(xCoordL(i)) !-(Az(i,j+1,k)-Az(i,j,k))/dy
!                 faceyData(MAG_FACE_VAR,i,j,k)= 0.0!(xCoordR(i)*Az(i+1,j,k)-xCoordL(i)*Az(i,j,k))/dx/xCoord(i)

              endif
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
                     (xCoordR(i)*facexData(MAG_FACE_VAR,i+1,j,  k  ) - xCoordL(i)*facexData(MAG_FACE_VAR,i,j,k))/dx/xCoord(i) &
                   + (faceyData(MAG_FACE_VAR,i,  j+1,k  ) - faceyData(MAG_FACE_VAR,i,j,k))/dy
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



