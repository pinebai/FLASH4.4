!!****if* source/Simulation/SimulationMain/unitTest/PFFT_PoissonFD/Simulation_init
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
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!
!! ARGUMENTS
!!
!!   none
!!
!!***

subroutine Simulation_init()

  use Simulation_data 
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none

  

  call RuntimeParameters_get('xmin', sim_xMin)
  call RuntimeParameters_get('xmax', sim_xMax)
  call RuntimeParameters_get('ymin', sim_yMin)
  call RuntimeParameters_get('ymax', sim_yMax)
  call RuntimeParameters_get('zmin', sim_zMin)
  call RuntimeParameters_get('zmax', sim_zMax)
  
end subroutine Simulation_init
