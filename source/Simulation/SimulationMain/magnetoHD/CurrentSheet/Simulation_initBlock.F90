!!****if* source/Simulation/SimulationMain/magnetoHD/CurrentSheet/Simulation_initBlock
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
!! PARAMETERS
!! 
!!
!!***

!!REORDER(4):solnData, face[xy]Data

subroutine Simulation_initBlock(blockId)

  use Simulation_data, ONLY : sim_gCell,sim_U0,sim_B0,sim_beta,&
                              sim_smallX,sim_smallP,sim_gamma,sim_killdivb,&
                              sim_xMin, sim_xMax, sim_yMin, sim_yMax

  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
                             Grid_getCellCoords,     &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr
  implicit none

#include "constants.h"
#include "Flash.h"

  ! compute the maximum length of a vector in each coordinate direction 
  ! (including guardcells)

  !!$ Arguments -----------------------
  integer, intent(in) :: blockId
  !!$ ---------------------------------

  integer :: i,j,k,n,istat,sizeX,sizeY,sizeZ
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real :: enerZone, ekinZone, eintZone
  real :: x1,x2,omega
  real, allocatable,dimension(:) :: xCoord,yCoord,zCoord
  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData

  
  ! dump some output to stdout listing the paramters
!!$   if (sim_meshMe == MASTER_PE) then
!!$1    format (1X, 1P, 4(A7, E13.7, :, 1X))
!!$2    format (1X, 1P, 2(A7, E13.7, 1X), A7, I13)
!!$  endif
  

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

  if (NDIM == 3) call Grid_getCellCoords&
                      (KAXIS, blockId, CENTER,sim_gCell, zCoord, sizeZ)
  if (NDIM >= 2) call Grid_getCellCoords&
                      (JAXIS, blockId, CENTER,sim_gCell, yCoord, sizeY)
  call Grid_getCellCoords(IAXIS, blockId, CENTER, sim_gCell, xCoord, sizeX)
  !------------------------------------------------------------------------------

  call Grid_getBlkPtr(blockID,solnData,CENTER)

#if NFACE_VARS > 0
  if (sim_killdivb) then
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
  endif
#endif


  x1 = (sim_xMax - sim_xMin)*0.25
  x2 = 3.*x1 + sim_xMin
  x1 = x1 + sim_xMin
  omega = 2*PI/(sim_yMax - sim_yMin)


  ! Loop over cells in the block.
  do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)


           ! Multiple species
           !solnData(SPECIES_BEGIN,i,j,k)=1.0e0-(NSPECIES-1)*sim_smallX
           do n=SPECIES_BEGIN,SPECIES_END
              solnData(n,i,j,k)=sim_smallX
           enddo

           ! Cell-centered values
           solnData(DENS_VAR,i,j,k)= 1.
           solnData(VELX_VAR,i,j,k)= sim_U0*sin(omega*yCoord(j))
           solnData(VELY_VAR,i,j,k)= 0.
           solnData(VELZ_VAR,i,j,k)= 0.
           solnData(PRES_VAR,i,j,k)=.5*sim_beta*sim_B0**2

           solnData(MAGX_VAR,i,j,k)= 0.
           solnData(MAGZ_VAR,i,j,k)= 0.

           if (xCoord(i) <= x1 ) then
              solnData(MAGY_VAR,i,j,k) =  sim_B0
              solnData(VECZ_VAR,i,j,k) = -solnData(MAGY_VAR,i,j,k)*xCoord(i) + x1
#if NFACE_VARS > 0
              if (sim_killdivb) then
                 facexData(MAG_FACE_VAR,i,j,k) = 0.
                 faceyData(MAG_FACE_VAR,i,j,k) = sim_B0
              endif
#endif
           elseif ((xCoord(i) >= x1 ).and.(xCoord(i) <= x2)) then
              solnData(MAGY_VAR,i,j,k) = -sim_B0
              solnData(VECZ_VAR,i,j,k) = -solnData(MAGY_VAR,i,j,k)*xCoord(i) - x1
#if NFACE_VARS > 0
              if (sim_killdivb) then
                 facexData(MAG_FACE_VAR,i,j,k) = 0.
                 faceyData(MAG_FACE_VAR,i,j,k) = -sim_B0
              endif
#endif
           else
              solnData(MAGY_VAR,i,j,k) =  sim_B0
              solnData(VECZ_VAR,i,j,k) = -solnData(MAGY_VAR,i,j,k)*xCoord(i) + 2.*x2-x1
#if NFACE_VARS > 0
              if (sim_killdivb) then
                 facexData(MAG_FACE_VAR,i,j,k) = 0.
                 faceyData(MAG_FACE_VAR,i,j,k) = sim_B0
              endif
#endif
           endif

           solnData(MAGP_VAR,i,j,k)=  .5*dot_product(solnData(MAGX_VAR:MAGZ_VAR,i,j,k),&
                                                     solnData(MAGX_VAR:MAGZ_VAR,i,j,k))
           solnData(DIVB_VAR,i,j,k)= 0.

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

#ifdef BETA_VAR
           solnData(BETA_VAR,i,j,k)=solnData(PRES_VAR,i,j,k)/solnData(MAGP_VAR,i,j,k)
#endif
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
  deallocate(yCoord)
  deallocate(zCoord)

end subroutine Simulation_initBlock










