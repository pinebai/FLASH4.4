!!****if* source/Simulation/SimulationMain/unitTest/Grid/UGReordered/Driver_evolveFlash
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
!!  This is the main global driver for simulations that are:
!!      Spatially refined, State form, strang split
!!
!! NOTES
!!
!!***

subroutine Driver_evolveFlash()

#include "constants.h"

  use Driver_data, ONLY: dr_globalMe 
  use Grid_interface, ONLY : Grid_unitTest

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

  print *, "Preparing to call Grid_unitTest"
  call Grid_unitTest(fileUnit,perfect)
  print *, "    Returned from Grid_unitTest"

  if (perfect) then
    write(fileUnit,'("all results conformed with expected values.")')
  endif

  close(fileUnit)

  return
  
end subroutine Driver_evolveFlash
