!!****if* source/Simulation/SimulationMain/unitTest/RungeKutta/2Dellipse/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!! SYNOPSIS
!!
!!  Simulation_init ()
!!
!! DESCRIPTION
!!
!!  Initializes the parameters for the Runge Kutta 2D ellipse unit test.
!!
!!***

subroutine Simulation_init ()

  use  Simulation_data
  use  RuntimeParameters_interface, ONLY: RuntimeParameters_get
  use  Logfile_interface,           ONLY: Logfile_stamp
  use  Driver_interface,            ONLY: Driver_abortFlash 
  use  sim_interface,               ONLY: sim_calculateInitialData

  implicit none

#include "constants.h"
#include "Flash.h"

  call RuntimeParameters_get ('sim_errorFraction',            sim_errorFraction          )
  call RuntimeParameters_get ('sim_RungeKuttaMethod',         sim_RungeKuttaMethod       )
  call RuntimeParameters_get ('sim_x0',                       sim_x0                     )
  call RuntimeParameters_get ('sim_y0',                       sim_y0                     )
  call RuntimeParameters_get ('sim_ellipseAspectRatio',       sim_ellipseAspectRatio     )
  call RuntimeParameters_get ('sim_stepSize',                 sim_stepSize               )
  call RuntimeParameters_get ('sim_numberOfEllipses',         sim_numberOfEllipses       )
!
!
!    ...Check parameters.
!
!
  if (sim_x0 /= 1.0 .or. sim_y0 /= 1.0) then
      call Driver_abortFlash ('[Simulation_init] ERROR: Particle position not at (1,1)!')
  end if

  if (sim_stepSize <= 0.0) then
      call Driver_abortFlash ('[Simulation_init] ERROR: Step size =< 0 !')
  end if
!
!
!    ...Set conversion factor(s).
!
!
  sim_deg2rad = acos (-1.0) / 180.0
!
!
!    ...Calculates the initital data.
!
!
  call sim_calculateInitialData ()
!
!
!    ...Ready!
!
!
  return
end subroutine Simulation_init
