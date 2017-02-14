!!****if* source/Simulation/SimulationMain/FlameChannel/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!! SYNOPSIS
!!
!!  call Simulation_init()
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   No arguments
!!
!! AUTOGENROBODOC
!!
!!
!!***

!  Initialize private data for the 3-stage flame test setup
!
! Dean Townsley 2008
!

subroutine Simulation_init()

  use Simulation_data
  use Flame_interface, ONLY : Flame_rhJump, Flame_getWidth, Flame_laminarSpeed
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use IO_interface, ONLY: IO_getScalar
  use Driver_data, ONLY: dr_restart
  use fl_effData, ONLY : fl_effDeltae, fl_eff_ye_u, fl_eff_sumy_u, &
     fl_eff_ye_b, fl_eff_sumy_b

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"

  real :: laminarWidth
  real :: ye_f, ye_a, yi_f, yi_a, qbar_f, qbar_a
  real, dimension(2) :: info
  character(len=MAX_STRING_LENGTH) :: xr_bcString
  integer :: xr_boundary_type

  integer :: i, ierr

  !--------------------------------------------------------
  !  initialize runtime parameters and some other constants
  !--------------------------------------------------------
  call RuntimeParameters_get( 'rho_ambient', sim_rhoAmbient)
  call RuntimeParameters_get( 't_ambient', sim_tAmbient)
  
  call RuntimeParameters_get( 'ignite', sim_ignite)
  call RuntimeParameters_get( 'frac_perturb', sim_fracPerturb)
  
  call RuntimeParameters_get( 'variableInflow', sim_variableInflow)
  call RuntimeParameters_get( 'inflowVortex', sim_inflowVortex)
  call RuntimeParameters_get( 'sigT', sim_sigT)
  call RuntimeParameters_get( 'sigP', sim_sigP)
  call RuntimeParameters_get( 'sigVx', sim_sigVx)
  call RuntimeParameters_get( 'sigVy', sim_sigVy)
  call RuntimeParameters_get( 'sigVz', sim_sigVz)

  ! sanitize sig's
  if ( sim_sigT < 0.0 ) then
     sim_sigT = 0.0
  else if ( sim_sigT > 1.0 ) then
     sim_sigT = 1.0
  endif

  if ( sim_sigP < 0.0 ) then
     sim_sigP = 0.0
  else if ( sim_sigP > 1.0 ) then
     sim_sigP = 1.0
  endif

  if ( sim_sigVx < 0.0 ) then
     sim_sigVx = 0.0
  else if ( sim_sigVx > 1.0 ) then
     sim_sigVx = 1.0
  endif

  if ( sim_sigVy < 0.0 ) then
     sim_sigVy = 0.0
  else if ( sim_sigVy > 1.0 ) then
     sim_sigVy = 1.0
  endif

  if ( sim_sigVz < 0.0 ) then
     sim_sigVz = 0.0
  else if ( sim_sigVz > 1.0 ) then
     sim_sigVz = 1.0
  endif

  call RuntimeParameters_get( 'yctr_vortex', sim_yctrVortex)
  call RuntimeParameters_get( 'vortexStrength', sim_vortexStrength)
  call RuntimeParameters_get( 'vortexSize', sim_vortexSize)
  call RuntimeParameters_get( 'xbegin_vortex', sim_xBeginVortex)
  call RuntimeParameters_get( 'xend_vortex', sim_xEndVortex)
  call RuntimeParameters_get( 'restart_vortex', sim_restartVortex)

  call RuntimeParameters_get('smooth_level', sim_smooth_level)
  call RuntimeParameters_get('vrms', sim_vrms)
  call RuntimeParameters_get('turbfield_filename',sim_turbfield_filename)
  call RuntimeParameters_get('turbfield_xmin',sim_turbfield_bbox(IAXIS,LOW))
  call RuntimeParameters_get('turbfield_xmax',sim_turbfield_bbox(IAXIS,HIGH))
  call RuntimeParameters_get('turbfield_ymin',sim_turbfield_bbox(JAXIS,LOW))
  call RuntimeParameters_get('turbfield_ymax',sim_turbfield_bbox(JAXIS,HIGH))
  call RuntimeParameters_get('turbfield_zmin',sim_turbfield_bbox(KAXIS,LOW))
  call RuntimeParameters_get('turbfield_zmax',sim_turbfield_bbox(KAXIS,HIGH))


  !  this is grid info, no accessor functions available
  call RuntimeParameters_get( 'xmin', sim_xmin)
  call RuntimeParameters_get( 'xmax', sim_xmax)
  call RuntimeParameters_get( 'ymin', sim_ymin)
  call RuntimeParameters_get( 'ymax', sim_ymax)
  call RuntimeParameters_get( 'zmin', sim_zmin)
  call RuntimeParameters_get( 'zmax', sim_zmax)

  sim_crossArea = (sim_ymax - sim_ymin) * (sim_zmax - sim_zmin)
 
  ! only need to get width of artificial flame once
  call Flame_getWidth(sim_laminarWidth)

  !--------------------------------------------------------
  !  find unburned and burned states
  !--------------------------------------------------------

  ! put this information in an eos datastructure and save qbar
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

  ! now determine properties of final NSE burned state
  call Flame_rhJump(sim_eosData_u, sim_eosData_b, fl_effDeltae, sim_flamespeed, MODE_DENS_TEMP)

  if ( .not. sim_ignite ) sim_eosData_b(:) = sim_eosData_u(:)

  if (dr_restart) then

     call IO_getScalar("last_burned_mass", sim_last_burned_mass)

  else

     sim_last_burned_mass = -1.0

  endif

  if (dr_restart .AND. sim_variableInflow) then

     call IO_getScalar("inflowVx", sim_inflowVx)

  else

     call RuntimeParameters_get("xr_boundary_type", xr_bcString)
     call RuntimeParameters_mapStrToInt(xr_bcString,xr_boundary_type)

     if (xr_boundary_type == USER_DEFINED) then ! assume it's inflow
        sim_inflowVx = -sim_flamespeed
     else ! we want the ash state to be zero
        !velx = sim_inflowVx - sim_flamespeed * ( alpha - 1.0 )
        !0.0 = sim_inflowVx - sim_flamespeed * (alpha - 1.0)
        sim_inflowVx = sim_flamespeed * &
           ( sim_eosData_u(EOS_DENS) / sim_eosData_b(EOS_DENS) - 1.0 )
     endif

  endif

end subroutine Simulation_init
