!!****if* source/Simulation/SimulationMain/unitTest/Opacity/Driver_initFlash
!!
!! NAME
!!  Driver_initFlash
!!
!! SYNOPSIS
!!
!!  Driver_initFlash ()
!!
!! DESCRIPTION
!!
!!  Stripped down version of the original Driver_initFlash.
!!  For comments on what is done please see the original routine.
!!
!! NOTES
!!
!!***

subroutine Driver_initFlash ()
  
  use Driver_data,                 ONLY : dr_globalComm,                 &
                                          dr_globalMe,                   &
                                          dr_globalNumProcs,             &
                                          dr_restart,                    &
                                          dr_elapsedWCTime,              &
                                          dr_initialWCTime,              &
                                          dr_particlesInitialized
  use Driver_interface,            ONLY : Driver_setupParallelEnv,       &
                                          Driver_init,                   &
                                          Driver_abortFlash
  use RuntimeParameters_interface, ONLY : RuntimeParameters_init
  use PhysicalConstants_interface, ONLY : PhysicalConstants_init
  use Logfile_interface,           ONLY : Logfile_init
  use Simulation_interface,        ONLY : Simulation_init
  use IO_interface,                ONLY : IO_init
  use Multispecies_interface,      ONLY : Multispecies_init
  use Grid_interface,              ONLY : Grid_init,                     &
                                          Grid_initDomain
  use RadTrans_interface,          ONLY : RadTrans_init
  use Opacity_interface,           ONLY : Opacity_init
  use Timers_interface,            ONLY : Timers_init,                   &
                                          Timers_start,                  &
                                          Timers_stop
  implicit none       
  
#include "constants.h"
#include "Flash.h"

  dr_elapsedWCTime = 0.0

  if (dr_globalNumProcs > 1) then
      call Driver_abortFlash ('[Driver_initFlash] ERROR: # of processors > 1')
  end if

  call dr_set_rlimits (dr_globalMe)

  call RuntimeParameters_init (dr_restart)
  write (*,*) ' Runtime parameters initialized'

  call Driver_setupParallelEnv ()
  write (*,*) ' Parallel environment set up'

  call Timers_init  (dr_initialWCTime)
  call Timers_start ("initialization")

  call PhysicalConstants_init ()
  write (*,*) ' Physical constants initialized'

  call Multispecies_init ()
  write (*,*) ' Multispecies initialized'

  call Logfile_init ()
  write (*,*) ' Logfile initialized'

  call Grid_init ()
  write (*,*) ' Grid initialized'

  call Driver_init ()
  write (*,*) ' Driver initialized'

  call Eos_init()
  write (*,*) ' Eos initialized'

  call IO_init ()
  write (*,*) ' IO initialized'

  call RadTrans_init ()
  write (*,*) ' RadTrans initialized'

  call Opacity_init ()
  write (*,*) ' Opacity initialized'

  call Grid_initDomain (dr_restart, dr_particlesInitialized)
  write (*,*) ' Grid domain initialized'

  call Simulation_init ()
  write (*,*) ' Simulation initialized'

  call Timers_stop ("initialization")

  return
end subroutine Driver_initFlash
