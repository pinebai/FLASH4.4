!!****if* source/Simulation/SimulationMain/unitTest/Pipeline/Driver_finalizeFlash
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
!!  Stripped down version from original for the Pipeline unit test.
!!
!!***
subroutine Driver_finalizeFlash ()

  use RuntimeParameters_interface, ONLY : RuntimeParameters_finalize
  use Grid_interface,              ONLY : Grid_finalize
  use Simulation_interface,        ONLY : Simulation_finalize
  use IO_interface,                ONLY : IO_finalize
  use Pipeline_interface,          ONLY : Pipeline_finalize
  use Timers_interface,            ONLY : Timers_finalize

  implicit none

#include "mpif.h"

  integer :: ierr
!
!
!     ...Proceed.
!
! 
  call RuntimeParameters_finalize ()
  call Grid_finalize              ()
  call Pipeline_finalize          ()
  call Simulation_finalize        ()
  call Timers_finalize            ()
  call IO_finalize                ()
  call MPI_Finalize               (ierr)

  return
end subroutine Driver_finalizeFlash
