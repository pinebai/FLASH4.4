!!****if* source/Simulation/SimulationMain/unitTest/Pipeline/Linear/sim_doAnalysis
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
!!  This routine analyzes the pipeline results. All items should have been transported
!!  (stored) to the last node with the highest rank. All other nodes should have no
!!  stored items.
!!
!! ARGUMENTS
!!
!!  perfect : the success indicator for the simulation
!!
!!***

subroutine sim_doAnalysis (perfect)

  use Simulation_Data,         ONLY : sim_baseName,             &
                                      sim_globalComm,           &
                                      sim_globalMe,             &
                                      sim_globalNumProcs,       &
                                      sim_lowestNumItemsOnProc, &
                                      sim_storedItems
                                      
  use Driver_interface,        ONLY : Driver_abortFlash

  implicit none

#include "Flash.h"
#include "constants.h"
 include "Flash_mpi.h"

  logical, intent (out) :: perfect

  integer :: highestRank
!
!
!     ...Act according to rank of node.
!
!
  write (*,*) ' processor ',sim_globalMe, ' has # of stored items = ',sim_storedItems

  highestRank = sim_globalNumProcs - 1

  if (sim_globalMe /= highestRank) then
      perfect = (sim_storedItems == 0)
  else
      perfect = (sim_storedItems == sim_lowestNumItemsOnProc)
  end if
!
!
!     ...Ready!
!
!
  return
end subroutine sim_doAnalysis
