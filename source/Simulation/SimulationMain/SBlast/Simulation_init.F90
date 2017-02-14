!!****if* source/Simulation/SimulationMain/SBlast/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  call Simulation_init()
!!
!!
!! DESCRIPTION
!!
!!  Initializes all the parameters needed for the Sedov-type blast
!!  problem
!!
!! ARGUMENTS
!!
!!  none
!!
!! PARAMETERS
!!
!!
!!***

subroutine Simulation_init()
  
  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, ONLY : Logfile_stamp
  implicit none

  
  
  call RuntimeParameters_get('sim_geo', sim_geo)
  call RuntimeParameters_get('sim_atmos1', sim_atmos1)
  call RuntimeParameters_get('sim_atmos2', sim_atmos2)

  call RuntimeParameters_get('sim_useE', sim_useE)
  call RuntimeParameters_get('sim_ibound', sim_ibound)

  call RuntimeParameters_get('sim_rhoIn', sim_rhoIn)
  call RuntimeParameters_get('sim_pIn', sim_pIn)
  call RuntimeParameters_get('sim_EIn', sim_EIn)
  call RuntimeParameters_get('sim_xcIn', sim_xcIn)
  call RuntimeParameters_get('sim_ycIn', sim_ycIn)
  call RuntimeParameters_get('sim_zcIn', sim_zcIn)
  call RuntimeParameters_get('sim_rIn', sim_rIn)

  call RuntimeParameters_get('sim_rho1', sim_rho1)
  call RuntimeParameters_get('sim_p1', sim_p1)
  call RuntimeParameters_get('sim_h1', sim_h1)
  call RuntimeParameters_get('sim_sh1', sim_sh1)

  call RuntimeParameters_get('sim_rho2', sim_rho2)
  call RuntimeParameters_get('sim_p2', sim_p2)
  call RuntimeParameters_get('sim_sh2', sim_sh2)

  call RuntimeParameters_get('xmin', xmin)
  call RuntimeParameters_get('xmax', xmax)
  call RuntimeParameters_get('ymin', ymin)
  call RuntimeParameters_get('ymax', ymax)
  call RuntimeParameters_get('zmin', zmin)
  call RuntimeParameters_get('zmax', zmax)

  call Logfile_stamp( "initializing Inhomogeneous Blast problem",  &
       "[Simulation_init]")
     
end subroutine Simulation_init







