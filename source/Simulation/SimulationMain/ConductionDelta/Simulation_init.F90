!!****if* source/Simulation/SimulationMain/ConductionDelta/Simulation_init
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
!!
!! DESCRIPTION
!!
!!  Initializes all the parameters needed for a particular simulation
!!
!!
!! ARGUMENTS
!!
!!  
!!
!! PARAMETERS
!!
!!  sim_xctr           Temperature peak center X-coordinate
!!  sim_yctr           Temperature peak center Y-coordinate
!!  sim_zctr           Temperature peak center Z-coordinate
!!***

subroutine Simulation_init()
  
  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  use Logfile_interface, ONLY : Logfile_stamp
  use Driver_interface, ONLY : Driver_getMype

  implicit none
#include "Flash.h"
#include "constants.h"

  

  real,external :: alogam
  real :: n ! sic
  real:: gasConstant,gammam1Inv
  integer :: iFault


  call RuntimeParameters_get('orientation', sim_orientation)
  call RuntimeParameters_get('rho_init', sim_rhoInit)
  call RuntimeParameters_get('toffset', sim_toffset)
  call RuntimeParameters_get('sim_Q', sim_Q)
  call RuntimeParameters_get('sim_tempBackground', sim_tempBackground)
  call RuntimeParameters_get('sim_xctr',sim_xCenter)
  call RuntimeParameters_get('sim_yctr',sim_yCenter)
  call RuntimeParameters_get('sim_zctr',sim_zCenter)
  call RuntimeParameters_get('smallp', sim_smallP)
  call RuntimeParameters_get('smallx', sim_smallX) 
  
  call RuntimeParameters_get('gamma', sim_gamma)
  call Driver_getMype(MESH_COMM, sim_meshMe)  

  call Logfile_stamp( "initializing Conduction problem",  &
       "[Simulation_init]")

  call RuntimeParameters_get("cond_TemperatureExponent", sim_condTemperatureExponent)
  call RuntimeParameters_get("iniCondTemperatureExponent", sim_iniCondTemperatureExponent)
  if (sim_iniCondTemperatureExponent==-999) then
     sim_iniCondTemperatureExponent = sim_condTemperatureExponent
  end if
  call RuntimeParameters_get("cond_K0", sim_alpha)
  call PhysicalConstants_get("ideal gas constant", gasConstant)
  gammam1Inv = 1.0/(sim_gamma-1.0)
  sim_alpha = sim_alpha / (sim_rhoInit * (gammam1Inv * gasConstant))
  n = sim_iniCondTemperatureExponent

  if (n .NE. 0) then
     sim_xi0 = (n+2)**(n+1) * 2.0**(1-n) / (n * (PI)**(0.5*n))
     sim_xi0 = sim_xi0 * exp(alogam(0.5+1.0/n,iFault)-alogam(1.0/n,iFault))
     sim_xi0 = sim_xi0**(1.0/(n+2))
     if (sim_meshMe == MASTER_PE) then
        print*,'Simulation_init: sim_xi0 is',sim_xi0
     end if


     sim_xfInitial = sim_xi0 * (sim_alpha * sim_Q**n * sim_toffset)**(1.0/(n+2))
     if (sim_meshMe == MASTER_PE) then
        print*,'Simulation_init: sim_xfInitial is',sim_xfInitial
     end if
  end if

end subroutine Simulation_init





