!!****if* source/Simulation/SimulationMain/unitTest/Opacity/Driver_evolveFlash
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
!!  Simple stripped down version for testing the opacity unit. 
!!
!! NOTES
!!
!!  The Driver unit uses a few unit scope variables that are
!!  accessible to all routines within the unit, but not to the
!!  routines outside the unit. These variables begin with "dr_"
!!  like, dr_globalMe or dr_dt, dr_beginStep, and are stored in fortran
!!  module Driver_data (in file Driver_data.F90. The other variables
!!  are local to the specific routine and do not have the prefix "dr_"
!!
!!
!!***

subroutine Driver_evolveFlash ()

  use Driver_data,         ONLY : dr_globalMe,                   &
                                  dr_nstep
  use Logfile_interface,   ONLY : Logfile_stamp,                 &
                                  Logfile_close
  use Timers_interface,    ONLY : Timers_start,                  &
                                  Timers_stop,                   &
                                  Timers_getSummary
  use Opacity_interface,   ONLY : Opacity_unitTest

  implicit none

# include "constants.h"
# include "Flash.h"

  character (len=20) :: fileName

  logical :: perfect
  integer :: temp,i

  integer, parameter     :: fileUnit = 2
  integer, dimension (4) :: prNum
!
!
!   ...Give a unique unitTest filename to each processor.
!
!  
  temp = dr_globalMe
  do i = 1,4
     prNum(i)= mod(temp,10)
     temp = temp/10
  end do
  filename = "unitTest_"//char(48+prNum(4))//char(48+prNum(3))//&
                                 char(48+prNum(2))//char(48+prNum(1))
  open  (fileUnit,file=fileName)
  write (fileUnit,'("P",I0)') dr_globalMe
!
!
!   ...Call the opacity unit test routine.
!
!  
  call Logfile_stamp ('Entering evolution routine' , '[Driver_evolveFlash]')
  call Timers_start  ("evolution")

  call Opacity_unitTest (fileUnit,perfect)

  if (perfect) then
     write (fileUnit,'("All results conformed with expected values.")')
  else
     write (fileUnit,'("Failure in Gravity unitTest at time")')
  end if
!
!
!   ...Final chores.
!
!  
  close (fileUnit)

  call Timers_stop   ("evolution")
  call Logfile_stamp ('Exiting evolution routine' , '[Driver_evolveFlash]')

  call Timers_getSummary (dr_nstep)

  call Logfile_close ()
!
!
!   ...Ready!
!
!  
  return
end subroutine Driver_evolveFlash



