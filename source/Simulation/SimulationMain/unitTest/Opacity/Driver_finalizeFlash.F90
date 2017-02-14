!!****if* source/Simulation/SimulationMain/unitTest/Opacity/Driver_finalizeFlash
!!
!! NAME
!!  Driver_finalizeFlash
!!
!! SYNOPSIS
!!  Driver_finalizeFlash ()
!!
!! DESCRIPTION
!!
!!  Stripped down version from original to test reading of opacity tables.
!!
!!***
subroutine Driver_finalizeFlash ()

  use RuntimeParameters_interface, ONLY : RuntimeParameters_finalize
  use Grid_interface,              ONLY : Grid_finalize
  use Simulation_interface,        ONLY : Simulation_finalize
  use IO_interface,                ONLY : IO_finalize
  use Multispecies_interface,      ONLY : Multispecies_finalize
  use RadTrans_interface,          ONLY : RadTrans_finalize
  use Opacity_interface,           ONLY : Opacity_finalize
  use Timers_interface,            ONLY : Timers_finalize

  implicit none

#include "mpif.h"

  integer :: ierr
 
  call RuntimeParameters_finalize ()

  call Grid_finalize ()
  write (*,*) ' Grid finalized'

  call Eos_finalize()             ! Equation of State
  write (*,*) ' Eos finalized'

  call Simulation_finalize ()
  write (*,*) ' Simulation finalized'

  call RadTrans_finalize ()
  write (*,*) ' RadTrans finalized'

  call Opacity_finalize ()
  write (*,*) ' Opacity finalized'

  call Multispecies_finalize ()
  write (*,*) ' Multispecies finalized'

  call Timers_finalize ()
  write (*,*) ' Timers finalized'

  call IO_finalize ()
  write (*,*) ' IO finalized'

  call MPI_Finalize (ierr)

  return
end subroutine Driver_finalizeFlash








