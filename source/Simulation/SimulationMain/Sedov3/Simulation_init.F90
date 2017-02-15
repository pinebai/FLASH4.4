!!****if* source/Simulation/SimulationMain/Sedov/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  Simulation_init()
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!  It calls RuntimeParameters_get rotuine for initialization.
!!  Initializes initial conditions for Sedov Spherical Explosion 
!!  problem.
!!
!! ARGUMENTS
!!
!!   
!!
!! PARAMETERS
!!
!!  sim_pAmbient       Initial ambient pressure
!!  sim_rhoAmbient     Initial ambient density
!!  sim_expEnergy      Explosion energy (distributed over 2^dimen central zones)
!!  sim_rInit          Radial position of inner edge of grid (for 1D )
!!  sim_xctr           Explosion center coordinates
!!  sim_yctr           Explosion center coordinates
!!  sim_zctr           Explosion center coordinates
!!  sim_nsubzones      Number of `sub-zones' in cells for applying 1d profile
!!
!!***

subroutine Simulation_init()

  use Simulation_data 
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface,           ONLY : Logfile_stamp

  implicit none

#include "constants.h"
#include "Flash.h"

  logical :: threadBlockListBuild, threadWithinBlockBuild  

  sim_pi = PI
  call RuntimeParameters_get('sim_pAmbient', sim_pAmbient)
  call RuntimeParameters_get('sim_rhoAmbient', sim_rhoAmbient)
  call RuntimeParameters_get('sim_expEnergy', sim_expEnergy)
  call RuntimeParameters_get('sim_rInit', sim_rInit)
  call RuntimeParameters_get('sim_nsubzones',sim_nSubZones)
  call RuntimeParameters_get('sim_xctr',sim_xCenter)
  call RuntimeParameters_get('sim_yctr',sim_yCenter)
  call RuntimeParameters_get('sim_zctr',sim_zCenter)
  call RuntimeParameters_get('gamma', sim_gamma)
  call RuntimeParameters_get('smallx', sim_smallX)
  call RuntimeParameters_get('smlrho', sim_smallRho)
  call RuntimeParameters_get('smallp', sim_smallP)
  call RuntimeParameters_get('smallT', sim_smallT)
  call RuntimeParameters_get('xmin',sim_xMin)
  call RuntimeParameters_get('ymin',sim_yMin)
  call RuntimeParameters_get('zmin',sim_zMin)
  call RuntimeParameters_get('xmax',sim_xMax)
  call RuntimeParameters_get('ymax',sim_yMax)
  call RuntimeParameters_get('zmax',sim_zMax)
  call RuntimeParameters_get('tinitial',sim_tInitial)

  !!!!!!!!!!!!!!!!!!!!!!!
  call RuntimeParameters_get('sim_expEnergy2', sim_expEnergy2)
  call RuntimeParameters_get('sim_xctr2',sim_xCenter2)
  call RuntimeParameters_get('sim_yctr2',sim_yCenter2)
  call RuntimeParameters_get('sim_zctr2',sim_zCenter2)
  !!!!!!!!!!!!!!!!!!!!!!!

  sim_iFuel = 1
  
  if (sim_nSubZones .le. 1) sim_nSubZones = 2
  
  sim_inSubZones = 1./real(sim_nSubZones)
  sim_inSubzm1   = 1./real(sim_nSubZones-1)
  sim_inszd      = sim_inSubZones**NDIM
  
  !
  !  Calculate the initial volume and interior pressure.
  !
  if (NDIM .eq. 1) then
     sim_vctr = 2. * sim_rInit
  elseif (NDIM .eq. 2) then
     sim_vctr = sim_pi * sim_rInit**2
  else
     sim_vctr = 4./3.*sim_pi*sim_rInit**3
  endif
  
  sim_pExp    = (sim_gamma-1.) * sim_expEnergy / sim_vctr
  sim_pExp2    = (sim_gamma-1.) * sim_expEnergy2 / sim_vctr



  call RuntimeParameters_get("threadBlockListBuild", threadBlockListBuild)
  call RuntimeParameters_get("threadHydroBlockList", sim_threadBlockList)

  call RuntimeParameters_get("threadWithinBlockBuild", threadWithinBlockBuild)
  call RuntimeParameters_get("threadHydroWithinBlock", sim_threadWithinBlock)

  if (sim_threadBlockList .and. .not. threadBlockListBuild) then
     call Logfile_stamp('WARNING! Turning off block list threading '//&
          'because FLASH is not built appropriately','[Simulation_init]')
     sim_threadBlockList = .false.
  end if
  if (sim_threadWithinBlock .and. .not. threadWithinBlockBuild) then
     call Logfile_stamp('WARNING! Turning off within block threading '//&
          'because FLASH is not built appropriately','[Simulation_init]')
     sim_threadWithinBlock = .false.
  end if
  
end subroutine Simulation_init
