!!****if* source/Simulation/SimulationMain/Flame1StageNoise/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!! SYNOPSIS
!!
!!  Simulation_init()
!!
!! DESCRIPTION
!!  Initialize private data for the 3-stage flame test setup
!!
!! ARGUMENTS
!!
!!
!!
!!***

subroutine Simulation_init()

  use Simulation_data
  use Flame_interface, ONLY : Flame_rhJump, Flame_getWidth, Flame_laminarSpeed
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use fl_effData, ONLY: fl_effDeltae, fl_eff_ye_u, fl_eff_sumy_u, fl_eff_ye_b, fl_eff_sumy_b

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"

  real :: laminarWidth

  !--------------------------------------------------------
  !  initialize runtime parameters and some other constants
  !--------------------------------------------------------
  call RuntimeParameters_get( 'rho_ambient', sim_rhoAmbient)
  call RuntimeParameters_get( 't_ambient', sim_tAmbient)
  
  call RuntimeParameters_get( 'ignite', sim_ignite)
  call RuntimeParameters_get( 'frac_perturb', sim_fracPerturb)
  
  call RuntimeParameters_get( 'pseudo_1d', sim_pseudo1d)
  call RuntimeParameters_get( 'xctr_perturb', sim_xctrPerturb)
  call RuntimeParameters_get( 'yctr_perturb', sim_yctrPerturb)
  call RuntimeParameters_get( 'zctr_perturb', sim_zctrPerturb)
  call RuntimeParameters_get( 'theta', sim_theta)
  
  !  this is grid info, no accessor functions available
  call RuntimeParameters_get( 'xmin', sim_xmin)
  call RuntimeParameters_get( 'xmax', sim_xmax)
  
  call RuntimeParameters_get( 'ymin', sim_ymin)
  call RuntimeParameters_get( 'ymax', sim_ymax)
  
  ! only need to get width of artificial flame once
  call Flame_getWidth(sim_laminarWidth)

  !--------------------------------------------------------
  !  find unburned and burned states
  !--------------------------------------------------------

  sim_eosData_u(EOS_DENS) = sim_rhoAmbient
  sim_eosData_u(EOS_TEMP) = sim_tAmbient
  sim_eosData_u(EOS_ABAR) = 1.e0 / fl_eff_sumy_u
  sim_eosData_u(EOS_ZBAR) = fl_eff_ye_u * sim_eosData_u(EOS_ABAR)

  sim_eosData_b(EOS_DENS) = sim_rhoAmbient
  sim_eosData_b(EOS_TEMP) = sim_tAmbient
  sim_eosData_b(EOS_ABAR) = 1.e0 / fl_eff_sumy_b
  sim_eosData_b(EOS_ZBAR) = fl_eff_ye_b * sim_eosData_b(EOS_ABAR)


  ! flamespeed should be constant
  call Flame_laminarSpeed(sim_eosData_u(EOS_DENS), sim_flamespeed)

  ! now determine praperties of final NSE burned state
  call Flame_rhJump(sim_eosData_u, sim_eosData_b, fl_effDeltae, sim_flamespeed, MODE_DENS_TEMP)


end subroutine Simulation_init
