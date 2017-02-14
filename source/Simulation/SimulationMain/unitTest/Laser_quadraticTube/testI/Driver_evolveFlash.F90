!!****if* source/Simulation/SimulationMain/unitTest/Laser_quadraticTube/testI/Driver_evolveFlash
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
!!  Simple stripped down version for testing the laser quadratic tube. 
!!
!! NOTES
!!
!!***

subroutine Driver_evolveFlash ()

  use Driver_data,       ONLY : dr_globalMe
  use Simulation_data,   ONLY : sim_printBlockVariables
  use Grid_interface,    ONLY : Grid_fillGuardCells
  use IO_interface,      ONLY : IO_outputFinal
  use Logfile_interface, ONLY : Logfile_close

  implicit none

# include "constants.h"
# include "Flash.h"

  character (len = 4                ) :: charProcessorID
  character (len = MAX_STRING_LENGTH) :: fileName

  logical :: perfect

  integer :: fileUnit
  integer :: ut_getFreeFileUnit
!
!
!   ...Open the (processor specific) indicator file. This is a file that will contain
!      the success (failure) status of the unit test.
!
!
  write (charProcessorID,'(I4.4)') dr_globalMe

  fileUnit = ut_getFreeFileUnit ()
  fileName = "unitTest_" // charProcessorID

  open (fileUnit, file = fileName)
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

  call sim_launchRays   ()
  call sim_doAnalysis   (perfect)
  call sim_printResults ()
!
!
!   ...Final chores. The exact phrase 'all results conformed with expected values.' must
!      be included to the indicator files to ensure recognition of a successful unit test
!      run by the 'flashTest/lib/flashModule.py' script.
!
!  
  if (perfect) then
      write (fileUnit,'(a)') 'SUCCESS all results conformed with expected values.'
  else
      write (fileUnit,'(a)') 'FAILURE'
  end if

  close (fileUnit)

  call IO_outputFinal ()

  call Logfile_close ()
!
!
!   ...Ready!
!
!  
  return
end subroutine Driver_evolveFlash



