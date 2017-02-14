!!****f* source/PhysicalConstants/PhysicalConstants_unitTest
!!
!! NAME
!!
!!  PhysicalConstants_unitTest
!!
!! SYNOPSIS
!!
!!  PhysicalConstants_unitTest(integer, intent(in)::fileUnit,
!!                             logical, intent(out)::perfect  )
!!
!! DESCRIPTION
!!
!!  This is a unitTest setup for testing the PhysicalConstants
!!      unit.  Normally called from Driver_evolveFlash within a Simulation unitTest
!!  See, for example, source/Simulation/SimulationMain/unitTest/PhysConst
!!
!! ARGUMENTS
!!      
!!     fileUnit - integer(in)   file for diagnostic output
!!     perfect  - logical(out)  TRUE if all tests are passed
!!
!! NOTES
!!
!!***

subroutine PhysicalConstants_unitTest(fileUnit,perfect)

  implicit none

  integer, intent(in)         :: fileUnit
  logical, intent(out)        :: perfect
  
  return 
  
end subroutine PhysicalConstants_unitTest
