!!****if* source/Simulation/SimulationMain/unitTest/Multispecies/Driver_evolveFlash
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
!!  This is the main global driver for testing the Multispecies
!!      unit
!!
!! NOTES
!!
!!***

subroutine Driver_evolveFlash()

  use Driver_data, ONLY: dr_globalMe
  use Multispecies_interface, ONLY : Multispecies_unitTest

#include "constants.h"

  implicit none

  integer                     ::  temp, prNum(4), i
  integer,parameter           ::  fileUnit = 2
  character(len=80)           ::  fileName
  logical ::  perfect = .true.

  !--------------------------------------------------------------

  ! stays true if no errors are found
  perfect = .true.

  ! Set up filename for parallel system
  temp = dr_globalMe
  do i = 1,4
     prNum(i)= mod(temp,10)
     temp = temp/10
  end do
  filename = "unitTest_"//char(48+prNum(4))//char(48+prNum(3))//&
       char(48+prNum(2))//char(48+prNum(1))

  open(fileUnit,file=fileName)

  if (dr_globalMe.EQ.MASTER_PE) print *,'Calling Multispecies_unitTest'
  call Multispecies_unitTest(fileUnit,perfect)
  if (dr_globalMe.EQ.MASTER_PE) print *,'Completed Multispecies_unitTest'


  if (perfect) then
     write(fileUnit,'("all results conformed with expected values.")')
  else
     write(fileUnit,'("errors found in test")')
  endif

  close(fileUnit)
  return 

end subroutine Driver_evolveFlash
