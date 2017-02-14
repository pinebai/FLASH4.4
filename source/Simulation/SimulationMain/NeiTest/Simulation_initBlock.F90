!!****if* source/Simulation/SimulationMain/NeiTest/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!! 
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer :: blockId, 
!!                       
!!
!!
!! DESCRIPTION
!!  
!!  Initializes a single block for the Neitest problem.
!!  This problem initializes a temperature step with plasma flow.
!!
!! ARGUMENTS
!!
!!  blockId -        The number of the block to initialize
!!  
!!
!!***
subroutine Simulation_initBlock(blockID)

#include "constants.h"
#include "Flash.h"
#include "Eos.h"

  use Simulation_speciesData, ONLY :  sim_specElementSymbol,sim_specElement,&
       sim_specSelected, sim_specNumElect, SPEC_NUM, SPEC_NIMAX, sim_specXfrac
  use Simulation_data, ONLY : sim_radius,sim_xstep, sim_smallx, sim_velInit,&
       sim_rhoAmbient, sim_tAmbient, sim_tPerturb, vel_init

  use Eos_interface, ONLY : Eos
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getCellCoords, &
       Grid_putRowData, Grid_putPointData
  use Ionize_interface, ONLY : Ionize_equil

  implicit none

  integer, intent(in) :: blockID
  

  integer i,j,k,n

  real,allocatable,dimension(:) :: xCenter,yCenter,zCenter
  integer,dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  integer :: iSize,jSize,kSize
  real, allocatable,dimension(:) :: rho, pres, temp, game, gamc, &
       vx, vy, vz, ener
  real, dimension(NSPECIES) :: massFrac
  logical :: gcell = .true.
  real :: xx, yy, zz
  real, dimension(EOS_NUM) :: eosData
  integer :: vecLen=1
  integer, dimension(1:MDIM) :: startingPos, cellIndex
  real :: temp_zone
  ! variables needed for the eos call
  real :: gamma, xalfa, xxni, xxne, xxnp


  real :: dist

  !=======================================================================


  integer :: id,nel,nion,inn 
  real             zel
  real,dimension(SPEC_NIMAX) :: delem

  real :: entropy, dst, dsd

  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)

  iSize = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
  jSize = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
  kSize = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1

  allocate(xCenter(iSize))
  allocate(yCenter(jSize))
  allocate(zCenter(kSize))
  allocate(rho(iSize))
  allocate(pres(iSize))
  allocate(temp(iSize)) 
  allocate(game(iSize))
  allocate(gamc(iSize))
  allocate(vx(iSize))
  allocate(vy(iSize))
  allocate(vz(iSize))
  allocate(ener(iSize))

  xCenter(:) = 0.0
  yCenter(:) = 0.0
  zCenter(:) = 0.0

  !Initialise vectors.  If we don't do this, we store uninitialised data in the 
  !guardcells.  This probably does not matter as the internal data is 
  !initialised, however, it stops uninitialised data run time warnings.
  rho(:) = 0.0; pres(:) = 0.0; temp(:) = 0.0; game(:) = 0.0; 
  gamc(:) = 0.0; vx(:) = 0.0; vy(:) = 0.0; vz(:) = 0.0; ener(:) = 0.0;
  

  if (NDIM == 3) call Grid_getCellCoords &
       (KAXIS, blockId, CENTER, gcell, zCenter, kSize)
  if (NDIM >= 2) call Grid_getCellCoords &
       (JAXIS, blockId, CENTER,gcell, yCenter, jSize)
  call Grid_getCellCoords(IAXIS, blockId, CENTER, gcell, xCenter, iSize)

  !-----------------------------------------------------------------------------
  ! loop over all the zones in the current block and set the temperature,
  ! density, and thermodynamics variables.
  !-----------------------------------------------------------------------------
  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     zz = zCenter(k)

     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        yy = yCenter(j)

        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           xx = xCenter(i)

           if (xx <= sim_xstep) then
              temp_zone = sim_tAmbient
           else
              temp_zone = sim_tPerturb
           endif

           !Copy in Electron & Hydrogen
           massFrac(1:2) = sim_specXfrac(1:2)

           id = 2
           do nel = 1,SPEC_NUM
              if (sim_specSelected(nel)) then
                 call Ionize_equil(temp_zone,nel,nion,delem)
                 nion = sim_specNumElect(nel)+1
                 do inn = 1,nion
                    id = id+1
                    massFrac(id) = delem(inn)*sim_specXfrac(id)
                 enddo
              endif
           enddo

           !
           !----------------------------------------------------------
           ! get the pressure and internal energy corresponding 
           ! to the ambient density
           ! and perturbed temperature
           !-----------------------------------------------------------

           eosData(:) = 0.0
           eosData(EOS_DENS)=sim_rhoAmbient
           eosData(EOS_TEMP)=temp_zone

           call Eos(MODE_DENS_TEMP,vecLen,eosData,massFrac)

           rho(i) = sim_rhoAmbient
           temp(i)   = temp_zone


           vx(i) = vel_init
           vy(i) = 0.0
           vz(i) = 0.0

           pres(i) = eosData(EOS_PRES)
           ener(i) = eosData(EOS_EINT) + 0.5*(vx(i)**2 + vy(i)**2 + vz(i)**2)

           game(i) = pres(i)/(eosData(EOS_EINT)*rho(i)) + 1.0
           gamc(i) = eosData(EOS_GAMC)

           !-----------------------------------------------------------------------------
           ! finally, fill the solution array
           !-----------------------------------------------------------------------------

           cellIndex(IAXIS)=i; cellIndex(JAXIS)=j; cellIndex(KAXIS)=k
           do n = 1, NSPECIES
              call Grid_putPointData(blockID, CENTER, SPECIES_BEGIN+n-1, EXTERIOR, &
                   cellIndex, massFrac(n))
           enddo

        enddo

        startingPos(IAXIS)=blkLimitsGC(LOW,IAXIS); startingPos(JAXIS)=j; startingPos(KAXIS)=k


        call Grid_putRowData(blockID,CENTER,DENS_VAR,EXTERIOR,IAXIS,&
             startingPos,rho,iSize)
        call Grid_putRowData(blockID,CENTER,ENER_VAR,EXTERIOR,IAXIS,&
             startingPos,ener,iSize)
        call Grid_putRowData(blockID,CENTER,TEMP_VAR,EXTERIOR,IAXIS,&
             startingPos,temp,iSize)
        call Grid_putRowData(blockID,CENTER,PRES_VAR,EXTERIOR,IAXIS,&
             startingPos,pres,iSize)
        call Grid_putRowData(blockID,CENTER,VELX_VAR,EXTERIOR,IAXIS,&
             startingPos,vx,iSize)
        call Grid_putRowData(blockID,CENTER,VELY_VAR,EXTERIOR,IAXIS,&
             startingPos,vy,iSize)
        call Grid_putRowData(blockID,CENTER,VELZ_VAR,EXTERIOR,IAXIS,&
             startingPos,vz,iSize)
        call Grid_putRowData(blockID,CENTER,GAME_VAR,EXTERIOR,IAXIS,&
             startingPos,game,iSize)
        call Grid_putRowData(blockID,CENTER,GAMC_VAR,EXTERIOR,IAXIS,&
             startingPos,gamc,iSize)

     enddo
  enddo

  deallocate(xCenter)
  deallocate(yCenter)
  deallocate(zCenter)
  deallocate(rho)
  deallocate(pres)
  deallocate(temp) 
  deallocate(game)
  deallocate(gamc)
  deallocate(vx)
  deallocate(vy)
  deallocate(vz)
  deallocate(ener)
  return
end subroutine Simulation_initBlock
