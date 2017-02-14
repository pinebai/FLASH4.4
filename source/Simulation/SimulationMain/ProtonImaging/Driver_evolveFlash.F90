!!****if* source/Simulation/SimulationMain/ProtonImaging/Driver_evolveFlash
!!
!! NAME
!!
!!  Driver_evolveFlash
!!
!! SYNOPSIS
!!
!!  Driver_evolveFlash ()
!!
!! DESCRIPTION
!!
!!  Simple stripped down version to prepare for calling the Proton Imaging unit. 
!!
!! NOTES
!!
!!***

subroutine Driver_evolveFlash ()

  use Driver_data,       ONLY : dr_nstep
  use Simulation_data,   ONLY : sim_printBlockVariables
  use Grid_interface,    ONLY : Grid_fillGuardCells
  use IO_interface,      ONLY : IO_output,IO_outputFinal
  use Timers_interface,  ONLY : Timers_getSummary
  use Logfile_interface, ONLY : Logfile_close

  implicit none

  logical endrun

# include "constants.h"
# include "Flash.h"
!
!
!     ...Fills all guard cells in the 'unk' array in all directions.
!
!
  call Grid_fillGuardCells (CENTER, ALLDIR)
!
!
!   ...Call the simulation routines.
!
!  
  if (sim_printBlockVariables) then
      call sim_printBlockData ()
  end if

  call sim_doProtonImaging ()
!
!
!   ...Final chores.
!
!  
!  call IO_output(0.0, 0.0, dr_nstep+1, 0, endrun, PLOTFILE_AND_PARTICLEFILE)
!  call IO_outputFinal ()

  call Timers_getSummary (dr_nstep)
  call Logfile_close ()
!
!
!   ...Ready!
!
!  
  return
end subroutine Driver_evolveFlash
