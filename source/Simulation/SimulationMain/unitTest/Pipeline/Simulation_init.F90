!!****if* source/Simulation/SimulationMain/unitTest/Pipeline/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!! SYNOPSIS
!!
!!  Simulation_init ()
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data for the pipeline
!!  unit test.
!!
!!***

subroutine Simulation_init()

  use Simulation_data

  use Driver_interface,            ONLY : Driver_abortFlash,  &
                                          Driver_getComm,     &
                                          Driver_getMype,     &
                                          Driver_getNumProcs

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none

#include "constants.h"
#include "Flash.h"

  integer :: item
  integer :: n, d
  integer :: sim_calculateMaxItems
!
!
!     ...Get the needed data.
!
!
  call Driver_getComm        (GLOBAL_COMM,                  sim_globalComm           )
  call Driver_getMype        (GLOBAL_COMM,                  sim_globalMe             )
  call Driver_getNumProcs    (GLOBAL_COMM,                  sim_globalNumProcs       )

  call RuntimeParameters_get ("basenm",                     sim_baseName             )
  call RuntimeParameters_get ("sim_itemSize",               sim_itemSize             )
  call RuntimeParameters_get ("sim_channelSize",            sim_channelSize          )
  call RuntimeParameters_get ("sim_lowestNumItemsOnProc",   sim_lowestNumItemsOnProc )
  call RuntimeParameters_get ("sim_maxItemsPipeline",       sim_maxItemsPipeline     )
!
!
!       ...Catch bad data.
!
!
  if (sim_itemSize < 1) then
      call Driver_abortFlash ('[Simulation_init] ERROR: Size of items must be > 0')
  end if

  if (sim_lowestNumItemsOnProc < 1) then
      call Driver_abortFlash ('[Simulation_init] ERROR: Lowest number of items/proc must be > 0')
  end if

  if (sim_channelSize < 1) then
      call Driver_abortFlash ('[Simulation_init] ERROR: Pipeline channel size must be > 0')
  end if
!
!
!       ...Do some needed chores.
!
!
  sim_baseName = adjustl (sim_baseName)    ! to get ready to use 'trim'
!
!
!       ...Calculate the depth of the Pascal triangle from the number of processors.
!          Determine the maximum number of items/proc (this would be at the beginning
!          of the Pascal triangle).
!
!
  n = sim_globalNumProcs
  d = int ((-3.0 + sqrt (8.0 * n + 1.0)) / 2)
  sim_depthPascalTriangle = d
  sim_cubeSide = int (real (n) ** (1.0/3.0))
  sim_maxItems = sim_calculateMaxItems ()
!
!
!       ...Allocate and set the items which will be sent through the pipeline.
!          Initially all items are located on node 0 and none are stored.
!
!
  allocate (sim_items (1:sim_itemSize, 1:sim_maxItems))

  if (sim_globalMe == 0) then

      do item = 1,sim_maxItems
         sim_items (1:sim_itemSize, item) = real (item)
      end do

      sim_numItems = sim_maxItems
  else
      sim_numItems = 0
  end if

  sim_storedItems = 0
!
!
!     ...Ready!
!
!
  return
end subroutine Simulation_init
