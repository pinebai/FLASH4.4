!!****if* source/Simulation/SimulationMain/Jeans/Simulation_init
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
!!  It calls RuntimeParameters_get rotuine for initialization.
!!  Initializes initial conditions for Jeans instability 
!!  problem.
!!
!! ARGUMENTS
!!
!!  none
!!
!! PARAMETERS
!!
!!  sim_p0       Initial ambient pressure
!!  sim_rho0     Initial ambient density
!!  sim_lambdaX  Perturbation X-wavelength
!!  sim_lambdaY  Perturbation Y-wavelength
!!  sim_lambdaZ  Perturbation Z-wavelength
!!  sim_A        Perturbation amplitude
!!
!!***

subroutine Simulation_init()

  use Simulation_data 
  use Driver_interface,            ONLY : Driver_getMype
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get

  implicit none

#include "constants.h"
#include "Flash.h"

  

  call Driver_getMype(MESH_COMM, sim_meshMe)

  sim_pi = PI
  call PhysicalConstants_get('Newton', sim_Newton)

  call RuntimeParameters_get('p0', sim_p0)
  call RuntimeParameters_get('rho0', sim_rho0)
  call RuntimeParameters_get('lambdax', sim_lambdaX)
  call RuntimeParameters_get('lambday', sim_lambdaY)
  call RuntimeParameters_get('lambdaz', sim_lambdaZ)
  call RuntimeParameters_get('amplitude', sim_A)
  call RuntimeParameters_get('gamma', sim_gamma)
  call RuntimeParameters_get('smalle', sim_smallE)
  call RuntimeParameters_get('delta_ref', sim_deltaRef)
  call RuntimeParameters_get('delta_deref', sim_deltaDeRef)
  call RuntimeParameters_get('reference_density', sim_refDensity)

  !
  !  Calculate the initial wavenumbers, sound speed, etc.
  !

  sim_kX = 2. * sim_pi / sim_lambdaX
  sim_kY = 2. * sim_pi / sim_lambdaY
  sim_kZ = 2. * sim_pi / sim_lambdaZ
  sim_kMag = sqrt(sim_kX**2 + sim_kY**2 + sim_kZ**2)
  sim_c0 = sqrt(sim_gamma*sim_p0/sim_rho0)
  sim_kJ = sqrt(4.*sim_pi*sim_Newton*sim_rho0)/sim_c0
  sim_oscFreq = (sim_c0*sim_kMag)**2 - 4.*sim_pi*sim_Newton*sim_rho0
  sim_velA = sim_A*sqrt(abs(sim_oscFreq)) / sim_kMag

  if (sim_meshMe == 0) then
     write (*,*)
     write (*,*) 'flash:  initializing for jeans problem.'
     write (*,*)
     write (*,*) 'kx   = ', sim_kX, 'ky   = ', sim_kY, 'kz   = ', sim_kZ
     write (*,*)
     write (*,*) 'k       = ', sim_kMag
     write (*,*) 'kJ      = ', sim_kJ
     write (*,*) 'c0      = ', sim_c0
     write (*,*) 'omega^2 = ', sim_oscFreq
     write (*,*) 'velamp  = ', sim_velA
     write (*,*)
  endif

  return

end subroutine Simulation_init
