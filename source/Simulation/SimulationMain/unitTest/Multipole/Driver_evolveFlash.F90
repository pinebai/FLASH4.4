!!****if* source/Simulation/SimulationMain/unitTest/Multipole/Driver_evolveFlash
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
!!  Simple stripped down version for testing the grid multipole solver. 
!!
!!***

subroutine Driver_evolveFlash ()

  use Driver_data,         ONLY : dr_globalMe,        &
                                  dr_nstep
  use Logfile_interface,   ONLY : Logfile_stamp,      &
                                  Logfile_close
  use Timers_interface,    ONLY : Timers_start,       &
                                  Timers_stop,        &
                                  Timers_getSummary
  use IO_interface,        ONLY : IO_output,          &
                                  IO_outputFinal

  implicit none

#include "constants.h"
#include "Flash.h"

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
!   ...Solve the poisson equation for all leaf blocks and analyze the obtained potential.
!
!  
  call Logfile_stamp ('Entering evolution routine' , '[Driver_evolveFlash]')
  call Timers_start  ("evolution")

  call sim_solveGridPoisson ()
  call sim_analyzePotential (perfect)
!
!
!   ...Final chores. The exact phrase 'all results conformed with expected values.' must
!      be included to the indicator files to ensure recognition of a successful unit test
!      run by the 'flashTest/lib/flashMudule.py' script.
!
!  
  if (perfect) then
      write (fileUnit,'(a)') 'SUCCESS all results conformed with expected values.'
  else
      write (fileUnit,'(a)') 'FAILURE'
  end if
!
!
!   ...Final chores.
!
!  
  close (fileUnit)

  call Timers_stop       ("evolution")
  call Logfile_stamp     ('Exiting evolution routine' , '[Driver_evolveFlash]')
  call IO_outputFinal    ()
  call Timers_getSummary (dr_nstep)
  call Logfile_close     ()


  return
  
end subroutine Driver_evolveFlash
