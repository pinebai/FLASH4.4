!!****if* source/Simulation/SimulationMain/magnetoHD/Rotor/Simulation_initBlock
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
!!  Reference: Balsara and Spicer, JCP, 149:270--292, 1999
!!             Balsara, The Astrophys Suppl Series, 151:148--184, 2004
!!
!! 
!! ARGUMENTS
!!
!!  blockID -          the number of the block to update
!!  
!!
!! 
!!
!!***

!!REORDER(4): solnData, face[xy]Data

subroutine Simulation_initBlock(blockID)

  use Simulation_data
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
                             Grid_getCellCoords,     &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr

  implicit none

#include "constants.h"
#include "Flash.h"

  !!$ Input Arguments --------------
  integer, intent(in) :: blockID
  !!$ ------------------------------

  integer :: i, j, k, n, istat, sizeX, sizeY, sizeZ
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real :: enerZone, ekinZone, eintZone
  real :: dx, dy, taper, radius, r0, v0
  real, allocatable, dimension(:) :: xCoord,yCoord,zCoord
  real, pointer, dimension(:,:,:,:) :: solnData, facexData, faceyData
  real :: perturbZ


  
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


  if (NDIM == 3) call Grid_getCellCoords(KAXIS,blockID,CENTER,sim_gCell,zCoord,sizeZ)
  if (NDIM >= 2) call Grid_getCellCoords(JAXIS,blockID,CENTER,sim_gCell,yCoord,sizeY)
  call Grid_getCellCoords(IAXIS,blockID,CENTER,sim_gCell,xCoord,sizeX)
  !------------------------------------------------------------------------------

  call Grid_getBlkPtr(blockID,solnData,CENTER)

#if NFACE_VARS > 0
  !if (sim_killdivb) then
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
  !endif
#endif

  ! Initialization for parameters
  r0=0.1
  v0=2.0

  ! Loop over cells in the block.
  do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)

           if (NSPECIES > 0) then
              ! Multiple species
              !solnData(SPECIES_BEGIN,i,j,k)=1.0e0-(NSPECIES-1)*sim_smallX
              do n=SPECIES_BEGIN,SPECIES_END
                 solnData(n,i,j,k)=sim_smallX
              enddo
           endif


           radius = sqrt((xCoord(i)-sim_xCtr)**2 + (yCoord(j)-sim_yCtr)**2)
           taper=(sim_Radius-radius)/(sim_Radius-r0)

           if (radius < r0) then
              solnData(DENS_VAR,i,j,k)=10.
              solnData(VELX_VAR,i,j,k)=-v0*(yCoord(j)-sim_yCtr)/.1
              solnData(VELY_VAR,i,j,k)= v0*(xCoord(i)-sim_xCtr)/.1
           elseif ((radius > r0).and.(radius < sim_Radius)) then
              solnData(DENS_VAR,i,j,k)=1.+9.*taper
              solnData(VELX_VAR,i,j,k)=-taper*v0*(yCoord(j)-sim_yCtr)/radius
              solnData(VELY_VAR,i,j,k)= taper*v0*(xCoord(i)-sim_xCtr)/radius
           else
              solnData(DENS_VAR,i,j,k)=1.
              solnData(VELX_VAR,i,j,k)=0.
              solnData(VELY_VAR,i,j,k)=0.
           endif

           solnData(VELZ_VAR,i,j,k) = 0.

           if (NDIM == 3) then
              ! velocity perturbations for 3D
              solnData(VELX_VAR,i,j,k)=solnData(VELX_VAR,i,j,k)*(1.+sim_perturbZ*sin(2.*PI*zCoord(k)))
              solnData(VELY_VAR,i,j,k)=solnData(VELY_VAR,i,j,k)*(1.+sim_perturbZ*sin(2.*PI*zCoord(k)))
              solnData(VELZ_VAR,i,j,k)=sim_perturbZ*sin(2.*PI*zCoord(k))
           endif

           solnData(DIVB_VAR,i,j,k) = 0.
           solnData(PRES_VAR,i,j,k) = 1.



           !Magnetic fields
           solnData(MAGX_VAR,i,j,k) = 5./sqrt(4.*PI)
           solnData(MAGY_VAR,i,j,k) = 0.
           solnData(MAGZ_VAR,i,j,k) = 0.
           solnData(MAGP_VAR,i,j,k)=  .5*dot_product(solnData(MAGX_VAR:MAGZ_VAR,i,j,k),&
                                                     solnData(MAGX_VAR:MAGZ_VAR,i,j,k))


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

           ! Cell face-centered variables for StaggeredMesh scheme
#if NFACE_VARS > 0
           if (sim_killdivb) then
              facexData(MAG_FACE_VAR,i,j,k)= 5./sqrt(4.*PI)
              faceyData(MAG_FACE_VAR,i,j,k)= 0.
           endif
#endif
        enddo
     enddo
  enddo


#if NFACE_VARS > 0
           if (sim_killdivb) then

              do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
                 i=blkLimitsGC(HIGH,IAXIS)+1
                 do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
                    facexData(MAG_FACE_VAR,i,j,k)= 5./sqrt(4.*PI)
                 enddo

                 j=blkLimitsGC(HIGH,JAXIS)+1
                 do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
                    faceyData(MAG_FACE_VAR,i,j,k)=0.
                 enddo
              enddo

           endif
#endif


  ! Release pointer
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)

#if NFACE_VARS > 0
  !if (sim_killdivb) then
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
  !endif
#endif


  deallocate(xCoord)
  deallocate(yCoord)
  deallocate(zCoord)


end subroutine Simulation_initBlock
