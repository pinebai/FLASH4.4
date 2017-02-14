!!****if* source/Simulation/SimulationMain/Chemistry_Test/Simulation_initBlock
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
!! DESCRIPTION
!!
!!
!!   Takes the one-d profile defined in Simulation_init and applys it
!!  This version sets up a spherical cluster in hydrostatic equilibrium
!!  in a 2D/3D cartesian grid. 
!!  The density and temperature profiles are give analytically.
!!
!!
!! ARGUMENTS
!!
!!  blockID -           the number of the block to update
!!
!!
!! PARAMETERS
!!
!!
!!***

subroutine Simulation_initblock (blockID)



!!***used modules from FLASH.***

  use Simulation_data
  use Logfile_interface, ONLY : Logfile_stamp
  use Grid_interface, ONLY : Grid_getCellCoords, Grid_getBlkPtr, &
       Grid_releaseBlkPtr, Grid_getBlkIndexLimits, Grid_putPointData
  
  use Eos_interface, ONLY : Eos 
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Multispecies.h"
#include "Eos.h"
  
  ! HERE are the arguments
  integer, intent(in) :: blockID
  
  
  
  !!***block number and (cell) counters***
  integer  ::  i, j, k, n
  
  
  !!** this is where we read and write the data
  real, pointer, dimension(:,:,:,:)  :: solnData
  
  !!***vectors that store cell dimensions **
  real, allocatable, dimension(:) :: x, y, z
  
  !!***coordinates of grid cells***
  real :: xx,  yy, zz
  
  !!***variables that you set as you go 
  real   ::   vx, vy, vz, p, rho, metals
  real   ::   e, ek, gam, T
  
  !!*** sizes in each of the three direction
  integer :: sizeX,sizeY,sizeZ
  integer :: istat
  
  !! This says that we are grabing the guard cells when we grab the block
  logical :: gcell = .true.
  
  
  
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  
  !This is a part dedicated to the PrimordialChemistry. This is followed
  !from the Cellular problem and the unitTest problem
  real :: abar, zbar
  real, dimension(SPECIES_BEGIN:SPECIES_END) :: massFraction
  real :: temp_zone, rho_zone, vel_zone
  real :: gamma, smallx
  real :: ptot, eint, etot
  real, dimension(EOS_NUM) :: eosData
  
  
  
  !  print *,'Simulation InitBlock'
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  call Grid_getBlkPtr(blockID,solnData)
  
  !!***get coordinates of this block***
  
  sizeX = blkLimitsGC(HIGH,IAXIS)
  sizeY = blkLimitsGC(HIGH,JAXIS)
  sizeZ = blkLimitsGC(HIGH,KAXIS)
  allocate(x(sizeX),stat=istat)
  allocate(y(sizeY),stat=istat)
  allocate(z(sizeZ),stat=istat)
  x = 0.0
  y = 0.0
  z = 0.0
  
  
  if (NDIM==3) call Grid_getCellCoords(KAXIS,blockID,CENTER,gcell,z,sizeZ)
  call Grid_getCellCoords(JAXIS,blockID,CENTER,gcell,y,sizey)
  call Grid_getCellCoords(IAXIS,blockID,CENTER,gcell,x,sizex)
  
  !Setup initial composition of all species
  !Question, Cellular had some small fraction for each species, can this not be zero?
  smallx = 1.0e-20;  !Trying something small
  massFraction(:) = smallx
  
  if (H_SPEC > 0) massFraction(H_SPEC) = max(sim_xH,smallx)
  if (HP_SPEC > 0) massFraction(HP_SPEC) = max(sim_xHP, smallx)
  if (HM_SPEC > 0) massFraction(HM_SPEC) = max(sim_xHM, smallx)
  if (D_SPEC > 0) massFraction(D_SPEC) = max(sim_xD, smallx)
  if (DP_SPEC > 0) massFraction(DP_SPEC) = max(sim_xDP, smallx)
  if (DM_SPEC > 0) massFraction(DM_SPEC) = max(sim_xDM, smallx)
  if (HE_SPEC > 0) massFraction(HE_SPEC) = max(sim_xHe, smallx)
  if (HEPP_SPEC > 0) massFraction(HEPP_SPEC) = max(sim_xHePP, smallx)
  if (HEP_SPEC > 0) massFraction(HEP_SPEC) = max(sim_xHeP, smallx)
  if (H2P_SPEC > 0) massFraction(H2P_SPEC) = max(sim_xH2P, smallx)
  if (H2_SPEC > 0) massFraction(H2_SPEC) = max(sim_xH2, smallx)
  if (HDP_SPEC > 0) massFraction(HDP_SPEC) = max(sim_xHDP, smallx)
  if (HD_SPEC > 0) massFraction(HD_SPEC) = max(sim_xHD, smallx)
  if (D2_SPEC > 0) massFraction(D2_SPEC) = max(sim_xD2, smallx)
  if (D2P_SPEC > 0) massFraction(D2P_SPEC) = max(sim_xD2P, smallx)
  if (ELEC_SPEC > 0) massFraction(ELEC_SPEC) = max(sim_xElec, smallx)
  
  
  !***************************************************
  !!***Now loop over all cells in the current block***
  !  
  print *,'Starting the InitBlock Loop'
  zz=0.
  do k = 1, sizeZ
     if(NDIM==3) zz = z(k)
     
     do j = 1, sizeY
        yy = y(j)
        
        do i = 1, sizeX
           xx = x(i)
           
           rho= sim_c_den 
           T= sim_c_temp 
           p=rho*T*sim_gasConst !!This will stay the same for pressure
           
           !Getting ready to find Gamma and some other Eos stuff
           eosData(EOS_TEMP) = T
           eosData(EOS_DENS) = rho
           eosData(EOS_PRES) = p
           
           call Eos(MODE_DENS_TEMP,1,eosData,massFraction)
           
           
           gamma = eosData(EOS_GAMC)
           sim_gamma = eosData(EOS_GAMC)
           ptot = eosData(EOS_PRES)
           eint = eosData(EOS_EINT)
           metals = 0.0
           
           vx   = 0.
           vy   = 0.
           vz   = 0.
           ek   = 0.
           
           
           solnData(DENS_VAR,i,j,k) = eosData(EOS_DENS)
           solnData(PRES_VAR,i,j,k) = eosData(EOS_PRES)
           solnData(VELX_VAR,i,j,k) = vx
           solnData(VELY_VAR,i,j,k) = vy
           solnData(VELZ_VAR,i,j,k) = vz
           solnData(GAME_VAR,i,j,k) = eosData(EOS_GAMC)
           solnData(GAMC_VAR,i,j,k) = eosData(EOS_GAMC)
           solnData(EINT_VAR,i,j,k) = eosData(EOS_EINT)
           solnData(ENER_VAR,i,j,k) = eosData(EOS_EINT) + ek
           solnData(TEMP_VAR,i,j,k) = eosData(EOS_TEMP)
           solnData(METL_MSCALAR,i,j,k) = sim_meta
           
           do n=SPECIES_BEGIN,SPECIES_END
              solnData(n, i,j,k) = massFraction(n)
           enddo
           
        enddo
     enddo
  enddo
  
  
  !print *, 'HERE'
  !  call Eos_wrapped(MODE_DENS_PRES,blkLimits,blockID)
  !  call Eos_wrapped(MODE_DENS_EI,blkLimits,blockID)
  !  call Eos_wrapped(MODE_DENS_TEMP,blkLimits,blockID)
  call Grid_releaseBlkPtr(blockID,solnData)
  
  !  print *, 'End of Simulation_initBlock.F90'
  
  deallocate(x)
  deallocate(y)
  deallocate(z)
  
  return
  
end subroutine Simulation_initBlock
