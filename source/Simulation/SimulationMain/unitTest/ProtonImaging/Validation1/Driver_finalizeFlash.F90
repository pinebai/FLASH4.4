!!****if* source/Simulation/SimulationMain/unitTest/ProtonImaging/Validation1/Driver_finalizeFlash
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
!!  Stripped down version from original for the Proton Imaging unit test.
!!
!!***
subroutine Driver_finalizeFlash ()

  use RuntimeParameters_interface, ONLY : RuntimeParameters_finalize
  use Grid_interface,              ONLY : Grid_finalize
  use Simulation_interface,        ONLY : Simulation_finalize
  use IO_interface,                ONLY : IO_finalize
  use ProtonImaging_interface,     ONLY : ProtonImaging_finalize
  use Timers_interface,            ONLY : Timers_finalize

  implicit none

#include "mpif.h"

  integer :: ierr
!
!
!     ...The 'Simulation_finalize' routine has to be called after
!        the 'ProtonImaging_finalize' routine, because it eventually
!        needs to read the proton detector files, which must be closed
!        at that moment.
!
! 
  call RuntimeParameters_finalize ()
  call Grid_finalize              ()
  call ProtonImaging_finalize     ()
  call Simulation_finalize        ()
  call Timers_finalize            ()
  call IO_finalize                ()
  call MPI_Finalize               (ierr)

  return
end subroutine Driver_finalizeFlash
