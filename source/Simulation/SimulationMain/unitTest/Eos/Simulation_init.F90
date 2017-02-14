!!****if* source/Simulation/SimulationMain/unitTest/Eos/Simulation_init
!!
!! NAME
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!  Simulation_init( )
!!
!!
!! DESCRIPTION
!!  Initializes all the parameters needed the Eos unit test
!!
!! ARGUMENTS
!!
!!  
!!
!!***

subroutine Simulation_init()
  
  use Simulation_data, ONLY : sim_xmin,sim_xmax,sim_ymin,sim_ymax,&
                              sim_zmin,sim_zmax,sim_smallx,&
                              sim_densMin,sim_tempMin,sim_xnMin,sim_presMin,&
                              sim_densMax, sim_tempMax, sim_xnMax, sim_presMax, &
                              sim_initialMass
  use Driver_interface, ONLY : Driver_abortFlash
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none

# include "Flash.h"

  

  integer :: lrefine_max, factor, nblockx,nblocky,nblockz

  call eos_initTest
  
  call RuntimeParameters_get( 'xmin', sim_xmin)
  call RuntimeParameters_get( 'xmax', sim_xmax)
  
  call RuntimeParameters_get( 'ymin', sim_ymin)
  call RuntimeParameters_get( 'ymax', sim_ymax)
  
  call RuntimeParameters_get( 'zmin', sim_zmin)
  call RuntimeParameters_get( 'zmax', sim_zmax)
  
  call RuntimeParameters_get( 'smallx', sim_smallx)

! sim_initialMass must be less than NSPECIES
  call RuntimeParameters_get( 'sim_initialMass', sim_initialMass)
  if (sim_initialMass .GT. NSPECIES)                        &
       call Driver_abortFlash('[Simulation_init] sim_initialMass must be less than NSPECIES')


  call RuntimeParameters_get( 'sim_densMin', sim_densMin)
  call RuntimeParameters_get( 'sim_densMax', sim_densMax)
  
  call RuntimeParameters_get( 'sim_tempMin', sim_tempMin)
  call RuntimeParameters_get( 'sim_tempMax', sim_tempMax)

  call RuntimeParameters_get( 'sim_presMin', sim_presMin)
  call RuntimeParameters_get( 'sim_presMax', sim_presMax)
  
  call RuntimeParameters_get( 'sim_xnMin', sim_xnMin)
  call RuntimeParameters_get( 'sim_xnMax', sim_xnMax)

end subroutine Simulation_init
