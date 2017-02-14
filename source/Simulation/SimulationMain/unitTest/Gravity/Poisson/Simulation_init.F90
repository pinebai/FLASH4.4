!!****if* source/Simulation/SimulationMain/unitTest/Gravity/Poisson/Simulation_init
!!
!! NAME
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!  Simulation_init( )
!!
!! DESCRIPTION   
!!   Initialize all the runtime parameters needed for the Particle unitTest
!!
!! ARGUMENTS
!!
!!   
!!
!! PARAMETERS
!!
!!   sim_rho_amb    Gas Density:  Entire domain receives this ambient parameter
!!   sim_p_amb      Gas Pressure: Entire domain receives this ambient parameter
!!   sim_vx_amb     Gas x-velocity:  Dominant flow velocity throughout domain 
!!   sim_vx_multiplier   Half of the domain in y has x-velocity multiplied by this value
!!   sim_seed   Random number seed -- NOT USED please ignore
!!   sim_vx_pert   Scales [-1,1] random number in x direction: set to zero for uniform flow
!!   sim_vy_pert   Scales [-1,1] random number in y direction: set to zero for uniform flow
!!   sim_vz_pert   Scales [-1,1] random number in z direction: set to zero for uniform flow
!!   sim_subSample  Subsampling of initial density distribution
!!
!!***

subroutine Simulation_init()

  use Simulation_data
  use Logfile_interface, ONLY : Logfile_stamp
  use Driver_interface,  ONLY : Driver_getMype
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none

#include "constants.h"
#include "Flash.h"
  
  
  integer :: i, j, status
 

  call Driver_getMype(MESH_COMM, sim_meshMe)

  !-------------------------------------------------------------------------------
  !               Write a message to stdout describing the problem setup.
  !-------------------------------------------------------------------------------
  

  if(sim_meshMe == MASTER_PE) then

     call Logfile_stamp( & 
          &         "initializing for Huang & Greengard problem", 'Simulation_init')
     write (*,*) & 
          &         'FLASH:  initializing for Huang & Greengard problem.'
     write (*,*)
     write (*,*) 'using ', Npeak, ' density peaks:'
     write (*,*)
     write (*,3) 'Peak', 'x', 'y', 'z', 'sigma'
     do i = 1, Npeak
        write (*,4) i, sim_xctr(i), sim_yctr(i), sim_zctr(i), sim_sigma(i)
     enddo
     write (*,*)
1    format (1X, 1P, 4(A13, E12.6, :, 1X))
2    format (1X, 1P, A13, I12)
3    format (1X, A5, 2X, 4(A13, :, 2X))
4    format (1X, I5, 2X, 1P, 4(E13.6, :, 2X))
  endif
  

  
  call RuntimeParameters_get('smlrho',sim_smlrho)
  call RuntimeParameters_get('sim_subSample', sim_subSample)
  

  return
end subroutine Simulation_init
