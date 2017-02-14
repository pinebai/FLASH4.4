!!****if* source/Simulation/SimulationMain/unitTest/Gravity/BHTree-cylinder/Simulation_init
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
!! ARGUMENTS
!!
!!    none
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!  It calls the RuntimeParameters_get routine for initialization.
!!
!!
!!***

subroutine Simulation_init()
  
  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  use Driver_interface, ONLY : Driver_abortFlash, Driver_getMype, Driver_getComm

  implicit none
#include "constants.h"


  call Driver_getMype(GLOBAL_COMM, sim_MyPE)
  call Driver_getComm(GLOBAL_COMM, sim_Comm)

  call RuntimeParameters_get( 'smlrho', smlrho)
  call RuntimeParameters_get( 'smallp', smallp)
  call RuntimeParameters_get( 'smallX', smallX)
  call RuntimeParameters_get( 'gamma_1', sim_gamma_1)
  call RuntimeParameters_get( 'gamma_1', gamma_1)
  call RuntimeParameters_get( 'gamma_2', gamma_2)
  call RuntimeParameters_get( 'abar_1', abar_1)
  call RuntimeParameters_get( 'abar_2', abar_2)
  call RuntimeParameters_get( 'xmax', xmax)
  call RuntimeParameters_get( 'xmin', xmin)
  call RuntimeParameters_get( 'ymax', ymax)
  call RuntimeParameters_get( 'ymin', ymin)
  call RuntimeParameters_get( 'zmax', zmax)
  call RuntimeParameters_get( 'zmin', zmin)

  call RuntimeParameters_get( 'sim_dens_c', sim_dens_c)
  call RuntimeParameters_get( 'sim_press_a', sim_press_a)
  call RuntimeParameters_get( 'sim_temp_c', sim_temp_c)
  call RuntimeParameters_get( 'sim_temp_a', sim_temp_a)
!  call RuntimeParameters_get( 'li_lambda', li_lambda)
  call RuntimeParameters_get( 'sim_solutionErrorTolerance1', sim_solutionErrorTolerance1)
  call RuntimeParameters_get( 'sim_solutionErrorTolerance2', sim_solutionErrorTolerance2)
  sim_Lx = xmax - xmin
  sim_Ly = ymax - ymin
  sim_Lz = zmax - zmin
 
  call PhysicalConstants_get( 'Boltzmann', sim_boltz)
  call PhysicalConstants_get( 'proton mass', sim_mH)
  call PhysicalConstants_get( 'Newton', sim_newt)
  call PhysicalConstants_get( 'pi', sim_pi)

  sim_vecLen = 1
  sim_mode = MODE_DENS_PRES

end subroutine Simulation_init
