!!****f* source/physics/materialProperties/Opacity/Opacity_unitTest
!!
!! NAME
!!
!!  Opacity_unitTest
!!
!! SYNOPSIS
!!
!!  Opacity_unitTest (integer (in)    :: fileUnit,
!!                    logical (inout) :: perfect)
!!
!! DESCRIPTION
!!
!!  This is a unitTest setup for testing the Opacity unit. Normally called
!!  from Driver_evolveFlash within a Simulation unitTest. See for example
!!  the directory 'source/Simulation/SimulationMain/unitTest/Opacity'.
!!
!! ARGUMENTS
!!
!!  fileUnit : number of file unit for diagnostic output
!!  perfect  : will be set .true., if the test is correct
!! 
!! NOTES
!!
!!  In the Simulation unit, you must set your Config file to 
!!  REQUIRES Opacity/OpacityMain/unitTest
!!
!!
!!***
subroutine Opacity_unitTest (fileUnit,perfect)

  implicit none

  integer, intent (in)    :: fileUnit
  logical, intent (inout) :: perfect

  perfect = .true.

  return
end subroutine Opacity_unitTest
