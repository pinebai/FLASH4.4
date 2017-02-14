!!****f* source/Multispecies/Multispecies_unitTest
!!
!! NAME
!!
!!  Multispecies_unitTest
!!
!! SYNOPSIS
!!
!!  Multispecies_unitTest(integer, intent(in)::fileUnit,
!!                        logical, intent(inout)::perfect  )
!!
!! DESCRIPTION
!!
!!  This is a unitTest setup for testing the Multispecies
!!      unit.  Normally called from Driver_evolveFlash within a Simulation unitTest
!!  See, for example, source/Simulation/SimulationMain/unitTest/Multispecies
!!
!! ARGUMENTS
!!
!!    fileUnit  -- number of file unit for diagnostic output
!!    perfect   -- if .true., unitTest has returned correctly
!! 
!! NOTES
!!  In the Simulation unit, you must set your Config file to 
!!    REQUIRES Multispecies/MultispeciesMain/unitTest
!!
!!***

subroutine Multispecies_unitTest(fileUnit,perfect)
!------------------------------------------------------------------------
    implicit none
    
 
    integer, intent(in)         :: fileUnit
    logical, intent(inout)      :: perfect
 
end subroutine Multispecies_unitTest
