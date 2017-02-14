!!****if* source/Simulation/SimulationMain/unitTest/RungeKutta/3Dcircle/Driver_initFlash
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
!!  Stripped down version for testing single units.
!!
!! NOTES
!!
!!***

subroutine Driver_initFlash ()
  
  use Driver_data,                 ONLY : dr_elapsedWCTime,          &
                                          dr_globalComm,             &
                                          dr_globalMe,               &
                                          dr_globalNumProcs,         &
                                          dr_initialWCTime,          &
                                          dr_particlesInitialized,   &
                                          dr_restart
  use Driver_interface,            ONLY : Driver_init,               &
                                          Driver_initNumericalTools, &
                                          Driver_setupParallelEnv
  use RuntimeParameters_interface, ONLY : RuntimeParameters_init
  use Simulation_interface,        ONLY : Simulation_init
  use Logfile_interface,           ONLY : Logfile_init

  implicit none       
  
  call dr_set_rlimits            (dr_globalMe)
  call RuntimeParameters_init    (dr_restart)
  call Driver_setupParallelEnv   ()
  call Timers_init               (dr_initialWCTime)
  call Logfile_init              ()
  call Driver_init               ()
  call Driver_initNumericalTools ()
  call Simulation_init           ()

  return
end subroutine Driver_initFlash
