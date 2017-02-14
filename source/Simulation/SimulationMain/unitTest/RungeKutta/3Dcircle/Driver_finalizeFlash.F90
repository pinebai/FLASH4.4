!!****if* source/Simulation/SimulationMain/unitTest/RungeKutta/3Dcircle/Driver_finalizeFlash
!!
!! NAME
!!
!!  Driver_finalizeFlash
!!
!! SYNOPSIS
!!
!!  Driver_finalizeFlash ()
!!
!! DESCRIPTION
!!
!!  Stripped down version from original for testing single units.
!!
!!***
subroutine Driver_finalizeFlash ()

  use RuntimeParameters_interface, ONLY : RuntimeParameters_finalize
  use Simulation_interface,        ONLY : Simulation_finalize
  use Timers_interface,            ONLY : Timers_finalize

  implicit none

#include "mpif.h"

  integer :: ierr
 
  call RuntimeParameters_finalize ()
  call Simulation_finalize ()
  call Timers_finalize ()
  call MPI_Finalize (ierr)

  return
end subroutine Driver_finalizeFlash
