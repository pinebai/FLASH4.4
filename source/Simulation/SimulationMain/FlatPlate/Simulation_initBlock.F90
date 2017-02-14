!!****if* source/Simulation/SimulationMain/FlatPlate/Simulation_initBlock
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
!!  Initializes fluid (initialized with -1.0) and 
!!  stationary rigid body data (initialized with 1.0)
!!
!!
!! ARGUMENTS
!!
!!  blockID -           the number of the block to update
!!***

!!REORDER(4): solnData

subroutine Simulation_initBlock(blockID)

  use Simulation_data
  use Grid_interface,  ONLY : Grid_getBlkIndexLimits, &
                              Grid_getBlkPtr, &
                              Grid_releaseBlkPtr, &
                              Grid_getCellCoords

  implicit none

#include "constants.h"
#include "Flash.h"

  integer,intent(IN) :: blockID
  real,pointer :: solnData(:,:,:,:)

  real :: rho_zone, velx_zone, vely_zone, velz_zone, pres_zone, &
       ener_zone, ekin_zone, eint_zone
  integer :: i, j, k, n, istat, sizeX, sizeY, sizeZ
  real,allocatable,dimension(:) :: xCoord,yCoord,zCoord
  real :: radius
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real, dimension(MDIM) :: del
  real :: naca_t,naca_c,naca_yt,naca_yc,naca_xc,naca_p,naca_m, naca_theta,xu,xl,yu,yl
  real :: soundSpeed, totalVel
!===============================================================================

  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)

  sizeX = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  sizeY = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
  sizeZ = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

  allocate(xCoord(sizeX),stat=istat)
  allocate(yCoord(sizeY),stat=istat)
  allocate(zCoord(sizeZ),stat=istat)

  xCoord = 0.0
  yCoord = 0.0
  zCoord = 0.0

  if (NDIM == 3) call Grid_getCellCoords(KAXIS,blockID,CENTER,sim_gCell,zCoord,sizeZ)
  if (NDIM >= 2) call Grid_getCellCoords(JAXIS,blockID,CENTER,sim_gCell,yCoord,sizeY)
  call Grid_getCellCoords(IAXIS,blockID,CENTER,sim_gCell,xCoord,sizeX)

  call Grid_getDeltas(blockID,del)



! In this problem the initial conditions are spatially uniform.
  
  
  rho_zone = sim_rhoAmbient
  pres_zone = sim_pAmbient


  soundSpeed = sqrt(sim_gamma*pres_zone/rho_zone)
  totalVel   = sim_Mach*soundSpeed



!!$  sim_windVelx = velx_zone
!!$  sim_windVely = vely_zone

  !! this is a 2D way of initializing x and y velocities,
  !! treat z velocity differently though.
  velx_zone = totalVel*cos(sim_xAngle)
  vely_zone = totalVel*sin(sim_xAngle)
  velz_zone = sim_windVelz
  
  ! Compute the gas energy and set the gamma-values needed for
  ! the equation of state.
  ekin_zone = 0.5 * (velx_zone**2 + vely_zone**2 + velz_zone**2)

  eint_zone = pres_zone / (sim_gamma-1.)
  eint_zone = eint_zone / rho_zone
  ener_zone = eint_zone + ekin_zone
  ener_zone = max(ener_zone, sim_smallP)
  

  call Grid_getBlkPtr(blockID, solnData, CENTER)
#if NSPECIES > 0
  solnData(SPECIES_BEGIN,:,:,:) =  1.0-(NSPECIES-1)*sim_smallX
  solnData(SPECIES_BEGIN+1:SPECIES_END,:,:,:) =     sim_smallX
#endif

  ! store the variables in the block's unk data
  solnData(DENS_VAR,:,:,:) = rho_zone
  solnData(PRES_VAR,:,:,:) = pres_zone
  solnData(ENER_VAR,:,:,:) = ener_zone
#ifdef EINT_VAR
  solnData(EINT_VAR,:,:,:) = eint_zone
#endif
  solnData(GAMC_VAR,:,:,:) = sim_gamma
  solnData(GAME_VAR,:,:,:) = sim_gamma


  solnData(VELX_VAR,:,:,:) = velx_zone
  solnData(VELY_VAR,:,:,:) = vely_zone
  solnData(VELZ_VAR,:,:,:) = velz_zone

  solnData(BDRY_VAR,:,:,:) = -1. !! -1 for fluid cells

  do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)

           if (sim_number == 1) then
              ! cylinder
              radius = (xCoord(i)-sim_xCtr)**2 + (yCoord(j)-sim_yCtr)**2
              radius = sqrt(radius)
              
              if (radius < sim_radius) then
                 solnData(DENS_VAR,i,j,k) = sim_rhoBulk
                 solnData(VELX_VAR,i,j,k) = 0.
                 solnData(VELY_VAR,i,j,k) = 0.
                 solnData(VELZ_VAR,i,j,k) = 0.
                 solnData(PRES_VAR,i,j,k) = pres_zone
                 solnData(BDRY_VAR,i,j,k) = 1. ! 1 for solid cells
              endif

           elseif (sim_number == 2) then

              ! NACA0015 airfoil
              if (xCoord(i) .ge. 0.0 .and. xCoord(i) .le. 1.0) then
                 naca_c=1
                 naca_t=0.15
                 naca_xc = xCoord(i)/naca_c
                 naca_yt=5.*naca_c*naca_t*( 0.2969*sqrt(naca_xc)&
                                           -0.1260*naca_xc&
                                           -0.3516*naca_xc**2&
                                           +0.2843*naca_xc**3&
                                           -0.1015*naca_xc**4 )

                 if (abs(yCoord(j)) .le. naca_yt) then
                    solnData(DENS_VAR,i,j,k) = sim_rhoBulk
                    solnData(VELX_VAR,i,j,k) = 0.
                    solnData(VELY_VAR,i,j,k) = 0.
                    solnData(VELZ_VAR,i,j,k) = 0.
                    solnData(PRES_VAR,i,j,k) = 2.*sim_pAmbient
                    solnData(BDRY_VAR,i,j,k) = 1.
                 endif
              endif

           elseif (sim_number == 3) then
              ! cambered NACA2412 airfoil
              naca_c=1
              naca_t=0.12
              naca_m=0.02
              naca_p=0.4
              naca_xc = xCoord(i)/naca_c

#if NDIM > 2
              if (abs(zCoord(k) - sim_zCtr) < 0.25) then
#endif
              if (xCoord(i) .ge. 0. .and. xCoord(i) .le. 1. ) then

                 naca_yt=5.*naca_c*naca_t*( 0.2969*sqrt(naca_xc)&
                                           -0.1260*naca_xc&
                                           -0.3516*naca_xc**2&
                                           +0.2843*naca_xc**3&
                                           -0.1015*naca_xc**4 )

                 if (xCoord(i) .ge. 0. .and. xCoord(i) .le. naca_p*naca_c) then
                    naca_yc= naca_m*(xCoord(i)/naca_p**2)*(2.*naca_p-naca_xc)

                    naca_theta=atan((naca_m/naca_p**2)*(2.*naca_p-naca_xc) &
                         - naca_m*xCoord(i)/(naca_p**2)/naca_c  )

                 elseif (xCoord(i) .gt. naca_p*naca_c .and. xCoord(i) .le. naca_c) then
                    naca_yc=naca_m*(naca_c-xCoord(i))/(1.-naca_p)**2*(1.+naca_xc-2.*naca_p)

                    naca_theta=atan((-naca_m/(1.-naca_p)**2)*(1.+naca_xc-2.*naca_p) &
                         + naca_m*(naca_c-xCoord(i))/naca_c/(1.-naca_p)**2)

                 endif

                 xu=xCoord(i)-naca_yt*sin(naca_theta)
                 xl=xCoord(i)+naca_yt*sin(naca_theta)
                 yu=naca_yc  +naca_yt*cos(naca_theta)
                 yl=naca_yc  -naca_yt*cos(naca_theta)

                 if (yCoord(j) .le. yu .and. yCoord(j) .ge. yl) then
                    solnData(DENS_VAR,i,j,k) = sim_rhoBulk
                    solnData(VELX_VAR,i,j,k) = 0.
                    solnData(VELY_VAR,i,j,k) = 0.
                    solnData(VELZ_VAR,i,j,k) = 0.
                    solnData(PRES_VAR,i,j,k) = 2.*sim_pAmbient
                    solnData(BDRY_VAR,i,j,k) = 1.
                 endif
              endif

#if NDIM > 2
              endif
#endif
              
           elseif (sim_number == 4) then
              !flat plate
              if (abs(xcoord(i)-sim_xCtr) < sim_radius  .and.  & !length of the plate
                  abs(yCoord(j)-sim_yCtr) < 2.*del(SWEEP_Y) ) then ! thickness of the plate
                 solnData(DENS_VAR,i,j,k) = sim_rhoBulk
                 solnData(VELX_VAR,i,j,k) = 0.
                 solnData(VELY_VAR,i,j,k) = 0.
                 solnData(VELZ_VAR,i,j,k) = 0.
                 solnData(PRES_VAR,i,j,k) = 2.*sim_pAmbient
                 solnData(BDRY_VAR,i,j,k) = 1.
              endif

           elseif (sim_number == 5) then
              ! solid wall in xCoord(i)>= 0.5
              if (xcoord(i) >= 0.5) then
                 solnData(DENS_VAR,i,j,k) = sim_rhoBulk
                 solnData(VELX_VAR,i,j,k) = 0.
                 solnData(VELY_VAR,i,j,k) = 0.
                 solnData(VELZ_VAR,i,j,k) = 0.
                 solnData(PRES_VAR,i,j,k) = 2.*sim_pAmbient
                 solnData(BDRY_VAR,i,j,k) = 1.
              endif


           elseif (sim_number == 6) then
              ! sphere
              radius = (xCoord(i)-sim_xCtr)**2 + (yCoord(j)-sim_yCtr)**2
#if NDIM > 2
              radius = radius + (zCoord(k)-sim_zCtr)**2
#endif
              radius = sqrt(radius)
              
              if (radius < sim_radius) then
                 solnData(DENS_VAR,i,j,k) = sim_rhoBulk
                 solnData(VELX_VAR,i,j,k) = 0.
                 solnData(VELY_VAR,i,j,k) = 0.
                 solnData(VELZ_VAR,i,j,k) = 0.
                 solnData(PRES_VAR,i,j,k) = pres_zone
                 solnData(BDRY_VAR,i,j,k) = 1. ! 1 for solid cells
              endif

           endif


        end do
     end do
  end do

  call Grid_releaseBlkPtr(blockID, solnData, CENTER)
 

  return
end subroutine Simulation_initBlock



