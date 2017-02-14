!!****if* source/Simulation/SimulationMain/unitTest/Pipeline/PascalTriangle2D/sim_doAnalysis
!!
!! NAME 
!!
!!  sim_doAnalysis
!!
!! SYNOPSIS
!!
!!  sim_doAnalysis (logical (out) :: perfect)
!!
!! DESCRIPTION
!!
!!  This routine analyzes the pipeline results.
!!
!! ARGUMENTS
!!
!!  perfect : the success indicator for the simulation
!!
!!***

subroutine sim_doAnalysis (perfect)

  use Simulation_Data,         ONLY : sim_baseName,            &
                                      sim_globalComm,          &
                                      sim_globalMe,            &
                                      sim_storedItems
                                      
  use Driver_interface,        ONLY : Driver_abortFlash

  implicit none

#include "Flash.h"
#include "constants.h"
 include "Flash_mpi.h"

  logical, intent (out) :: perfect
!
!
!     ...xxx
!
!
  write (*,*) ' processor ',sim_globalMe, ' has # of stored items = ',sim_storedItems
!
!
!     ...The 'perfect' indicator will be set initially to false as a default on all
!        processors.
!
!
  perfect = .false.
!
!
!     ...Ready!
!
!
  return
end subroutine sim_doAnalysis
