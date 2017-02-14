!!****if* source/Simulation/SimulationMain/StirTurb/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer, intent(IN)  :: blockid)
!!
!! DESCRIPTION
!!
!!  Initializes data (density, pressure, velocity, etc.) for
!!  a specified block.  This version sets up the Stirring-Turbulence
!!  problem.
!!
!!
!! ARGUMENTS
!!
!!   blockid : ID of block in current processor
!!
!!
!!
!!
!!***


subroutine Simulation_initBlock(blockId)

  use Simulation_data, ONLY: sim_rhoAmbient, sim_cAmbient, sim_gamma
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_getCellCoords, Grid_putPointData, Grid_putRowData
 
  implicit none 
#include "constants.h"
#include "Flash.h"
#include "Eos.h"

  integer, intent(IN) :: blockId
  real, dimension(NSPECIES) :: massFrac

#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_IHI_GC) :: xCoord
  real, dimension(GRID_JHI_GC) :: yCoord
  real, dimension(GRID_KHI_GC) :: zCoord
  real, dimension(GRID_IHI_GC) :: rho, pressure, energy, e, temperature
  real, dimension(GRID_IHI_GC) :: vx, vy, vz, gamE, gamC,scalar
  real, dimension(GRID_IHI_GC) :: accx, accy, accz
#else
  real,allocatable,dimension(:) :: xCoord,yCoord,zCoord
  real, dimension(:),allocatable :: rho, pressure, energy, e, temperature
  real, dimension(:),allocatable :: vx, vy, vz, gamE, gamC,scalar
  real, dimension(:),allocatable :: accx, accy, accz
#endif
  real :: yy, zz, p_ambient
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
  integer :: sizeX,sizeY,sizeZ,sizeIntX
  integer,dimension(MDIM) :: startingPos
  

  integer :: i, j, k, n,pass_scalar,istat
  logical :: gcell = .true.

  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
  sizeX = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  sizeIntX = blkLimits(HIGH,IAXIS)-blkLimits(LOW,IAXIS)+1
  sizeY = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
  sizeZ = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1


#ifndef FIXEDBLOCKSIZE
  allocate(xCoord(sizeX),stat=istat)
  allocate(yCoord(sizeY),stat=istat)
  allocate(zCoord(sizeZ),stat=istat)
  allocate(rho(sizeIntX),stat=istat)
  allocate(pressure(sizeIntX),stat=istat)
  allocate(energy(sizeIntX),stat=istat)
  allocate(e(sizeIntX),stat=istat)
  allocate(temperature(sizeIntX),stat=istat)
  allocate(vx(sizeIntX),stat=istat)
  allocate(vy(sizeIntX),stat=istat)
  allocate(vz(sizeIntX),stat=istat)
  allocate(gamE(sizeIntX),stat=istat)
  allocate(gamC(sizeIntX),stat=istat)
  allocate(scalar(sizeIntX),stat=istat)
  allocate(accx(sizeIntX),stat=istat)
  allocate(accy(sizeIntX),stat=istat)
  allocate(accz(sizeIntX),stat=istat)
#endif

  ! get the coordinate information for the current block from the database
  xCoord = 0.0
  yCoord = 0.0
  zCoord = 0.0

#if NSPECIES > 0
     massFrac = 0.0
     massFrac(1) = 1.0
#endif

  if (NDIM == 3) call Grid_getCellCoords&
                      (KAXIS, blockId, CENTER, gcell, zCoord, sizeZ)
  if (NDIM >= 2) call Grid_getCellCoords&
                      (JAXIS, blockId, CENTER, gcell, yCoord, sizeY)
  call Grid_getCellCoords(IAXIS, blockId, CENTER, gcell, xCoord, sizeX)

  pass_scalar = NPROP_VARS+NSPECIES+1
  do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
     zz = zCoord(k)

     do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
        yy = yCoord(j)

        p_ambient      = sim_rhoAmbient * sim_cAmbient**2 / sim_gamma
        rho        (:) = sim_rhoAmbient
        pressure   (:) = p_ambient
        energy     (:) = p_ambient / ( (sim_gamma - 1.)  * rho)
!        temperature(:) = sim_eosArr(EOS_TEMP)
        gamC       (:) = sim_gamma
        gamE           = sim_gamma
        vx         (:) = 0.
        vy         (:) = 0.
        vz         (:) = 0.

        accx       (:) = 0.
        accy       (:) = 0.
        accz       (:) = 0.
        scalar     (:) = zz
        do i = blkLimits(LOW,IAXIS), blkLimits(HIGH, IAXIS)

           startingPos(IAXIS) = i
           startingPos(JAXIS) = j
           startingPos(KAXIS) = k

           do n=1,NSPECIES
              call Grid_putPointData(blockId, CENTER, SPECIES_BEGIN+n-1, EXTERIOR, startingPos, massFrac(n))

           enddo  ! end loop over species
        enddo ! end loop over i
        startingPos(IAXIS)=blkLimits(LOW,IAXIS)
        e = energy + 0.5*(vx*vx + vy*vy + vz*vz)
 

        !this routine initializes a row at a time...

        call Grid_putRowData(blockId, CENTER, ACCX_VAR, EXTERIOR, IAXIS, startingPos, accx, sizeIntX)
        call Grid_putRowData(blockId, CENTER, ACCY_VAR, EXTERIOR, IAXIS, startingPos, accy, sizeIntX)
        call Grid_putRowData(blockId, CENTER, ACCZ_VAR, EXTERIOR, IAXIS, startingPos, accz, sizeIntX)

        call Grid_putRowData(blockId, CENTER, DENS_VAR, EXTERIOR, IAXIS, startingPos, rho, sizeIntX)
#ifdef EINT_VAR
        call Grid_putRowData(blockId, CENTER, EINT_VAR, EXTERIOR, IAXIS, startingPos, energy, sizeIntX)
#endif
        call Grid_putRowData(blockId, CENTER, ENER_VAR, EXTERIOR, IAXIS, startingPos, e, sizeIntX)
        call Grid_putRowData(blockId, CENTER, PRES_VAR, EXTERIOR, IAXIS, startingPos, pressure, sizeIntX)
        call Grid_putRowData(blockId, CENTER, TEMP_VAR, EXTERIOR, IAXIS, startingPos, temperature, sizeIntX)
        call Grid_putRowData(blockId, CENTER, GAME_VAR, EXTERIOR, IAXIS, startingPos, gamE, sizeIntX)
        call Grid_putRowData(blockId, CENTER, GAMC_VAR, EXTERIOR, IAXIS, startingPos, gamC, sizeIntX)
        call Grid_putRowData(blockId, CENTER, VELX_VAR, EXTERIOR, IAXIS, startingPos, vx, sizeIntX)
        call Grid_putRowData(blockId, CENTER, VELY_VAR, EXTERIOR, IAXIS, startingPos, vy, sizeIntX)
        call Grid_putRowData(blockId, CENTER, VELZ_VAR, EXTERIOR, IAXIS, startingPos, vz, sizeIntX)
        !call Grid_putRowData(blockId, CENTER, pass_scalar, EXTERIOR, IAXIS, startingPos, scalar, sizeIntX)

     enddo ! end loop over j

  enddo ! end loop over k
#ifndef FIXEDBLOCKSIZE
  deallocate(xCoord)
  deallocate(yCoord)
  deallocate(zCoord)
  deallocate(rho)
  deallocate(pressure)
  deallocate(energy)
  deallocate(e)
  deallocate(temperature)
  deallocate(vx)
  deallocate(vy)
  deallocate(vz)
  deallocate(gamE)
  deallocate(gamC)
  deallocate(scalar)
  deallocate(accx)
  deallocate(accy)
  deallocate(accz)
#endif
end subroutine Simulation_initBlock
