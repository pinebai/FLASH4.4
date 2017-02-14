!!****if* source/Simulation/SimulationMain/unitTest/RungeKutta/3Dcircle/Simulation_init
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
!!  Initializes the parameters for the Runge Kutta unit test.
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
  call RuntimeParameters_get ('sim_rx0',                      sim_rx0                    )
  call RuntimeParameters_get ('sim_ry0',                      sim_ry0                    )
  call RuntimeParameters_get ('sim_rz0',                      sim_rz0                    )
  call RuntimeParameters_get ('sim_speed',                    sim_speed                  )
  call RuntimeParameters_get ('sim_stepSize',                 sim_stepSize               )
  call RuntimeParameters_get ('sim_numberOfCircles',          sim_numberOfCircles        )
  call RuntimeParameters_get ('sim_numberOfRungeKuttaSteps',  sim_numberOfRungeKuttaSteps)
!
!
!    ...Check parameters.
!
!
  if (sim_rx0 == 0.0 .and. sim_ry0 == 0.0) then
      call Driver_abortFlash ('[Simulation_init] ERROR: Particle position on z-axis!')
  end if

  if (sim_speed <= 0.0) then
      call Driver_abortFlash ('[Simulation_init] ERROR: Particle speed =< 0 !')
  end if

  if (sim_stepSize <= 0.0) then
      call Driver_abortFlash ('[Simulation_init] ERROR: Step size =< 0 !')
  end if
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
