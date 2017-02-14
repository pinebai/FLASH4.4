!!****if* source/Simulation/SimulationMain/unitTest/Roots/x3Polynomials/Driver_evolveFlash
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
!!  Simple stripped down version for testing single units.
!!
!! NOTES
!!
!!***

subroutine Driver_evolveFlash ()

  use Driver_data,  ONLY : dr_globalMe

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
!   ...Do the x3 polynomials root test.
!
!  
  call sim_x3PolynomialsRootsTest (perfect)
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

  return
end subroutine Driver_evolveFlash
