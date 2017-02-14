!!****if* source/Simulation/SimulationMain/magnetoHD/OrszagTang/Simulation_initBlock
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
!!  Reference: Orszag and Tang, J. Fluid Mech., 90:129--143, 1979
!!
!! 
!! ARGUMENTS
!!
!!  blockID -          the number of the block to update
!!
!! 
!!
!!***

!!REORDER(4): solnData, face[xy]Data

subroutine Simulation_initBlock(blockID)

  use Simulation_data, ONLY : sim_gCell, sim_gamma,   &
                              sim_smallX, sim_smallP, &
                              sim_killdivb, sim_perturbation

  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
                             Grid_getCellCoords,     &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr

  implicit none

#include "constants.h"
#include "Flash.h"

  !!$ Arguments -----------------------
  integer, intent(in) :: blockID
  !!$ ---------------------------------

  integer :: i, j, k, n, istat, sizeX, sizeY, sizeZ

  real,allocatable,dimension(:) :: xCoord,yCoord,zCoord
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real :: enerZone, ekinZone, eintZone
  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData


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
  if (sim_killdivb) then
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
     facexData(:,:,:,:)=0.0
     faceyData(:,:,:,:)=0.0


     if (NDIM == 3) then
        call Grid_getBlkPtr(blockID,facezData,FACEZ)
        facezData(:,:,:,:) = 0.
     endif

  endif
#endif

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

           ! Cell-centered values
           if (NDIM == 2) then
              solnData(DENS_VAR,i,j,k)=  1.
              solnData(VELX_VAR,i,j,k)= -sin(2.*PI*yCoord(j))
              solnData(VELY_VAR,i,j,k)=  sin(2.*PI*xCoord(i))
              solnData(VELZ_VAR,i,j,k)=  0.
              solnData(PRES_VAR,i,j,k)=  1./sim_gamma

           elseif (NDIM == 3) then
              !! 3D setup by Helzel, Rossmanith, and Taetz, JCP 230, 2011
              solnData(DENS_VAR,i,j,k)=  sim_gamma**2
              solnData(VELX_VAR,i,j,k)= -(1.+sim_perturbation*sin(2.*PI*zCoord(k)))&
                   *sin(2.*PI*yCoord(j))
              solnData(VELY_VAR,i,j,k)=  (1.+sim_perturbation*sin(2.*PI*zCoord(k)))&
                   *sin(2.*PI*xCoord(i))
              solnData(VELZ_VAR,i,j,k)=  sim_perturbation*sin(2.*PI*zCoord(k))
              solnData(PRES_VAR,i,j,k)=  sim_gamma
           endif

           solnData(MAGX_VAR,i,j,k)= -sin(2.*PI*yCoord(j))
           solnData(MAGY_VAR,i,j,k)=  sin(4.*PI*xCoord(i))
           solnData(MAGZ_VAR,i,j,k)=  0.

           if (NDIM == 2) then
              solnData(MAGX_VAR:MAGZ_VAR,i,j,k) = solnData(MAGX_VAR:MAGZ_VAR,i,j,k)/sim_gamma
           endif


           solnData(MAGP_VAR,i,j,k)= .5*dot_product(solnData(MAGX_VAR:MAGZ_VAR,i,j,k),&
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
              facexData(MAG_FACE_VAR,i,j,k)=-sin(2.*PI*yCoord(j))
              faceyData(MAG_FACE_VAR,i,j,k)= sin(4.*PI*xCoord(i))
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
           facexData(MAG_FACE_VAR,i,j,k)=-sin(2.*PI*yCoord(j))
        enddo

        j=blkLimitsGC(HIGH,JAXIS)+1
        do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
           faceyData(MAG_FACE_VAR,i,j,k)= sin(4.*PI*xCoord(i))
        enddo
     enddo

  endif

  if (NDIM == 2) then
     facexData = facexData/sim_gamma
     faceyData = faceyData/sim_gamma
  elseif (NDIM == 3) then
     facezData(MAG_FACE_VAR,:,:,:) = 0.
  endif
#endif


  ! Release pointer

#if NFACE_VARS > 0
  if (sim_killdivb) then
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
  endif
#endif

  call Grid_releaseBlkPtr(blockID,solnData,CENTER)

  deallocate(xCoord)
  deallocate(yCoord)
  deallocate(zCoord)

end subroutine Simulation_initBlock






