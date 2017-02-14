!!****if* source/Simulation/SimulationMain/unitTest/PFFT_PoissonFD/Driver_evolveFlash
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
!!  This is the global driver specifically for the Grid unit test
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

subroutine Driver_evolveFlash()

#include "constants.h"

  use Driver_data, ONLY: dr_globalMe
  use Grid_interface, ONLY : Grid_unitTest
  use IO_interface, ONLY : IO_writeCheckpoint, IO_outputFinal

  implicit none

  logical ::  perfect = .true.

  character(len=20) :: fileName
  integer, parameter        :: fileUnit = 2
  integer,dimension(4) :: prNum
  integer :: temp,i


  ! stays true if no errors are found
  perfect = .true.

  temp = dr_globalMe

  do i = 1,4
     prNum(i)= mod(temp,10)
     temp = temp/10
  end do
  filename = "unitTest_"//char(48+prNum(4))//char(48+prNum(3))//&
                                 char(48+prNum(2))//char(48+prNum(1))

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

  close(fileUnit)

  ! Dumping results
  call IO_writeCheckpoint()
  call IO_outputFinal( )

  return
  
end subroutine Driver_evolveFlash
