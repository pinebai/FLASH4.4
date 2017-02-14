!!****if* source/Simulation/SimulationMain/unitTest/Poisson/Driver_evolveFlash
!!
!! NAME
!!
!!  Driver_evolveFlash
!!
!! SYNOPSIS
!!
!!  call Driver_evolveFlash()
!!
!! DESCRIPTION
!!
!!  This is the global driver specifically for a Poisson GridSolver  unit test.
!!  It opens a file for each of the processor that it is running
!!  on and call the Grid_unitTest with the file unit. The Grid_unitTest
!!  functions carries out the testing of the unit, and reports success or
!!  failure through a logical argument "perfect". Upon return from
!!  Grid_unitTest, the Driver_evolveFlash function makes appropriate notification
!!  in the file which the test suite can parse to determine if the test was
!!  successful.
!!
!!
!!***


#include "Flash.h"

subroutine Driver_evolveFlash()

  use Driver_data, ONLY: dr_globalMe, &
                         dr_nstep

  use Logfile_interface,ONLY : Logfile_stamp, Logfile_close
  use Timers_interface, ONLY : Timers_start, Timers_stop, &
                               Timers_getSummary

  use Grid_interface,    ONLY : Grid_unitTest !,        &
                                !Grid_getListOfBlocks

  use IO_interface,      ONLY : IO_outputFinal


  implicit none

#include "constants.h"

  integer :: blockCount
  integer :: blockList(MAXBLOCKS)


  logical ::  perfect

  character(len=20) :: fileName
  integer           :: fileUnit, ut_getFreeFileUnit
  integer,dimension(4) :: prNum
  integer :: temp,i


  ! stays true if no errors are found
  perfect = .true.

  call Logfile_stamp('Entering evolution loop' , '[Driver_evolveFlash]')
  call Timers_start("evolution")


  temp = dr_globalMe

  do i = 1,4
     prNum(i)= mod(temp,10)
     temp = temp/10
  end do
  filename = "unitTest_"//char(48+prNum(4))//char(48+prNum(3))//&
                                 char(48+prNum(2))//char(48+prNum(1))

  fileUnit = ut_getFreeFileUnit()
  open(fileUnit,file=fileName)
  write(fileUnit,'("P",I0)') dr_globalMe

  if (dr_globalMe .eq. 0) print *, "Preparing to call Grid_unitTest"
  call Grid_unitTest(fileUnit,perfect)
  if (dr_globalMe .eq. 0)   print *, "Returned from Grid_unitTest"


  if (perfect) then
    write(fileUnit,'("all results conformed with expected values.")')
  else
    write(fileUnit,'("test failed.")')
  endif

  call Timers_stop("evolution")
  call Logfile_stamp('Exiting evolution loop' , '[Driver_evolveFlash]')

  close(fileUnit)

  call IO_outputFinal()
  call Timers_getSummary(dr_nstep)
  call Logfile_stamp("FLASH run complete.", "LOGFILE_END")
  call Logfile_close()

  return

end subroutine Driver_evolveFlash
