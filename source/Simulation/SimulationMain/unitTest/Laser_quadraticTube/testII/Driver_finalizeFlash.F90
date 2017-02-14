!!****if* source/Simulation/SimulationMain/unitTest/Laser_quadraticTube/testII/Driver_finalizeFlash
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
!!  Stripped down version from original for the laser quadratic tube.
!!
!!***
subroutine Driver_finalizeFlash ()

  use RuntimeParameters_interface, ONLY : RuntimeParameters_finalize
  use Grid_interface,              ONLY : Grid_finalize
  use Simulation_interface,        ONLY : Simulation_finalize
  use IO_interface,                ONLY : IO_finalize
  use Eos_interface,               ONLY : Eos_finalize
  use EnergyDeposition_interface,  ONLY : EnergyDeposition_finalize
  use Timers_interface,            ONLY : Timers_finalize

  implicit none

#include "mpif.h"

  integer :: ierr
 
  call RuntimeParameters_finalize ()
  call Grid_finalize ()
  call Simulation_finalize ()
  call EnergyDeposition_finalize ()
  call Timers_finalize ()
  call IO_finalize ()
  call Eos_finalize ()
  call MPI_Finalize (ierr)

  return
end subroutine Driver_finalizeFlash
