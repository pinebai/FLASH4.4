!!****if* source/Simulation/SimulationMain/unitTest/PhysConst/Driver_evolveFlash
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
!!  This is the main global driver for testing the PhysicalConstants
!!      unit
!!  Tthe unitTest has no physical domain or time stepping,
!!      so the regular calls to grid, etc.
!!      are not needed.  Instead, it simply calls PhysicalConstants_unitTest
!!      which organizes all the initialization etc.
!!
!! NOTES
!!
!!***

subroutine Driver_evolveFlash()

  use PhysicalConstants_interface, ONLY : PhysicalConstants_unitTest

  implicit none

  integer                     ::  temp, prNum(4), i
  integer,parameter           ::  fileUnit = 2
  logical                     ::  perfect
  character(len=80)           ::  fileName
  
  !--------------------------------------------------------------
  
  ! stays true if no errors are found
  perfect = .true.
  
  do i = 1,4
     prNum(i)= mod(temp,10)
     temp = temp/10
  end do
  filename = "unitTest_"//char(48+prNum(4))//char(48+prNum(3))//&
                                 char(48+prNum(2))//char(48+prNum(1))

  open(fileUnit,file=fileName)
  
  call PhysicalConstants_unitTest(fileUnit,perfect)
    
  if (perfect) then
     write(fileUnit,'("all results conformed with expected values.")')
  else
     write(fileUnit,*)"FAILURE in unit test"
  endif
  
  close(fileUnit)
  return 
  
end subroutine Driver_evolveFlash
