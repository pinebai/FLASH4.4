!!****if* source/Simulation/SimulationMain/magnetoHD/BrioWu/Simulation_init
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
!!  none
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!  It calls RuntimeParameters_get routine for initialization.
!!  Initializes initial conditions for BrioWu problem.
!!
!!***

subroutine Simulation_init()

  use Simulation_data, ONLY : sim_xangle, sim_yangle, sim_smallx,  sim_smallP,   &
                              sim_posn,   sim_gamma,  sim_killdivb,&
                              sim_uRight, sim_uLeft,  sim_rhoRight,sim_rhoLeft,  &
                              sim_pRight, sim_pLeft,  sim_byLeft,  sim_byRight,  &
                              sim_bNormal,sim_gcell,  sim_xcos,    sim_ycos, sim_zcos,&
                              sim_vRight, sim_vLeft,  sim_wRight,  sim_wLeft,    &
                              sim_bzLeft, sim_bzRight,sim_xmin,    sim_xmax,     &
                              sim_ymin,   sim_ymax,   sim_meshMe


  use Driver_interface,  ONLY : Driver_getMype
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  use Logfile_interface, ONLY : Logfile_stamp

  implicit none

#include "constants.h"
#include "Flash.h"

  

  call Driver_getMype(MESH_COMM, sim_meshMe)

  call RuntimeParameters_get('xangle',    sim_xangle)
  call RuntimeParameters_get('yangle',    sim_yangle)
  call RuntimeParameters_get('smallx',    sim_smallx)
  call RuntimeParameters_get('smallp',    sim_smallP)
  call RuntimeParameters_get('posn',      sim_posn)
  call RuntimeParameters_get('gamma',     sim_gamma)
  call RuntimeParameters_get('rho_right', sim_rhoRight)
  call RuntimeParameters_get('rho_left',  sim_rhoLeft)
  call RuntimeParameters_get('u_right',   sim_uRight)
  call RuntimeParameters_get('u_left',    sim_uLeft)
  call RuntimeParameters_get('v_right',   sim_vRight)
  call RuntimeParameters_get('v_left',    sim_vLeft)
  call RuntimeParameters_get('w_right',   sim_wRight)
  call RuntimeParameters_get('w_left',    sim_wLeft)
  call RuntimeParameters_get('p_right',   sim_pRight)
  call RuntimeParameters_get('p_left',    sim_pLeft)
  call RuntimeParameters_get('by_left',   sim_byLeft)
  call RuntimeParameters_get('by_right',  sim_byRight)
  call RuntimeParameters_get('bz_left',   sim_bzLeft)
  call RuntimeParameters_get('bz_right',  sim_bzRight)
  call RuntimeParameters_get('b_normal',  sim_bNormal)
  call RuntimeParameters_get('killdivb',  sim_killdivb)
  call RuntimeParameters_get('xmin',      sim_xmin)
  call RuntimeParameters_get('xmax',      sim_xmax)
  call RuntimeParameters_get('ymin',      sim_ymin)
  call RuntimeParameters_get('ymax',      sim_ymax)

  if (NDIM == 1) then
     sim_killdivb = .false.
  endif

  sim_gcell = .true.

  ! Compute angle for the shock normal
  sim_xangle = sim_xangle * 0.0174532925        ! Convert to radians.
  sim_yangle = sim_yangle * 0.0174532925

  ! Compute direction cosines for the shock normal.       
  sim_xcos = cos(sim_xangle)
  if (NDIM .eq. 1) then
     sim_xcos = 1.
     sim_ycos = 1. !0.
     sim_zcos = 1. !0.
  elseif (NDIM .eq. 2) then
     sim_ycos = sqrt(1. - sim_xcos*sim_xcos)
     sim_zcos = 0.
  elseif (NDIM .eq. 3) then
     sim_ycos = cos(sim_yangle)
     sim_zcos = sqrt( max(0., 1. - sim_xcos*sim_xcos - sim_ycos*sim_ycos) )
  endif
     

end subroutine Simulation_init
