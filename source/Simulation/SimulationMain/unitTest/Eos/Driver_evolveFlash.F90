!!****if* source/Simulation/SimulationMain/unitTest/Eos/Driver_evolveFlash
!!
!! NAME
!!
!!  Driver_evolveFlash
!!
!! SYNOPSIS
!!
!!  Driver_evolveFlash()
!!
!! DESCRIPTION
!!
!!  This is the main for the Eos unit test
!!
!!***

subroutine Driver_evolveFlash()

#include "constants.h"

  use Driver_data, ONLY: dr_globalMe 
  use Eos_interface, ONLY : Eos_unitTest
  implicit none

  interface
     integer function ut_getFreeFileUnit()
     end function ut_getFreeFileUnit
  end interface

  logical :: perfect

  character(len=20) :: fileName
  integer           :: fileUnit
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

  fileUnit = ut_getFreeFileUnit()
  open(fileUnit,file=fileName)
  write(fileUnit,'("P",I0)') dr_globalMe

  print *, "Preparing to call Eos_unitTest"
  Call Eos_unitTest(fileUnit,perfect)
  print *, "    Returned from Eos_unitTest"

  if (perfect) then
    write(fileUnit,'("all results conformed with expected values.")')
  endif

  close(fileUnit)


  return
  
end subroutine Driver_evolveFlash
