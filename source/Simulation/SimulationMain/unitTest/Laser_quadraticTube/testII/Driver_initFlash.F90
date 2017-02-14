!!****if* source/Simulation/SimulationMain/unitTest/Laser_quadraticTube/testII/Driver_initFlash
!!
!! NAME
!!
!!  Driver_initFlash
!!
!! SYNOPSIS
!!
!!  Driver_initFlash ()
!!
!! DESCRIPTION
!!
!!  Stripped down version of the original for the laser quadratic tube.
!!  For comments on what is done please see the original routine.
!!
!! NOTES
!!
!!***

subroutine Driver_initFlash ()
  
  use Driver_data,                 ONLY : dr_elapsedWCTime,        &
                                          dr_globalComm,           &
                                          dr_globalMe,             &
                                          dr_globalNumProcs,       &
                                          dr_initialWCTime,        &
                                          dr_particlesInitialized, &
                                          dr_restart
  use Driver_interface,            ONLY : Driver_init,             &
                                          Driver_setupParallelEnv
  use RuntimeParameters_interface, ONLY : RuntimeParameters_init
  use PhysicalConstants_interface, ONLY : PhysicalConstants_init
  use Logfile_interface,           ONLY : Logfile_init
  use Simulation_interface,        ONLY : Simulation_init
  use IO_interface,                ONLY : IO_init,                 &
                                          IO_outputInitial
  use Eos_interface,               ONLY : Eos_init
  use EnergyDeposition_interface,  ONLY : EnergyDeposition_init
  use Grid_interface,              ONLY : Grid_init,               &
                                          Grid_initDomain
  use Timers_interface,            ONLY : Timers_init,             &
                                          Timers_start,            &
                                          Timers_stop
  implicit none       
  
#include "constants.h"
#include "Flash.h"

  dr_elapsedWCTime = 0.0

  call dr_set_rlimits          (dr_globalMe)
  call RuntimeParameters_init  (dr_restart)
  call Driver_setupParallelEnv ()
  call Timers_init             (dr_initialWCTime)
  call Timers_start            ("initialization")
  call PhysicalConstants_init  ()
  call Logfile_init            ()
  call Grid_init               ()
  call Driver_init             ()
  call IO_init                 ()
  call Eos_init                ()
  call EnergyDeposition_init   ()
  call Simulation_init         ()
  call Grid_initDomain         (dr_restart, dr_particlesInitialized)
  call IO_outputInitial        (1, 0.0)
  call Timers_stop             ("initialization")

  return
end subroutine Driver_initFlash
