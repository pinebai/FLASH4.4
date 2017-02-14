!!****if* source/Simulation/SimulationMain/unitTest/Pipeline/sim_doPipelineTest
!!
!!  NAME 
!!
!!   sim_doPipelineTest
!!
!!  SYNOPSIS
!!
!!   sim_doPipelineTest ()
!!
!!  DESCRIPTION
!!
!!   This routine transports all the items from node 0 through the pipeline.
!!
!! ARGUMENTS
!!
!!***

subroutine sim_doPipelineTest ()

  use Pipeline_interface,  ONLY : Pipeline_globalCheckStatus,    &
                                  Pipeline_globalCheckStructure, &
                                  Pipeline_localActivate,        &
                                  Pipeline_localCreate,          &
                                  Pipeline_localDeactivate,      &
                                  Pipeline_localDestroy,         &
                                  Pipeline_localFlush,           &
                                  Pipeline_localGetItems,        &
                                  Pipeline_localProgress

  use Simulation_data,     ONLY : sim_channelSize,      &
                                  sim_items,            &
                                  sim_itemSize,         &
                                  sim_maxItems,         &
                                  sim_maxItemsPipeline, &
                                  sim_numItems

  implicit none

#include "Flash.h"
#include "constants.h"

  logical :: activePipeline
  logical :: empty
  logical :: fullestChannelOnly
  logical :: localLoop
  logical :: moreItems

  integer :: numGetItems
!
!
!     ...Create/check/activate the pipeline.
!
!
  call Pipeline_localCreate (sim_itemSize,         &
                             sim_maxItemsPipeline, &
                             sim_channelSize,      &
                             'Test'                )

  call Pipeline_globalCheckStructure ()
  call Pipeline_localActivate        ()
!
!
!     ...Move items locally through the pipeline. If there are active items present
!        on the local processor, put them on the sending channels, while at the same
!        time receiving new items. Items are continued to be put on the sending
!        channels, until there are no more present. Then sending channels are
!        forcefully flushed, while transportation through the pipeline progresses.
!        Check global status of pipeline, if all items have been sent and all
!        sending channels are empty.
!
!
  activePipeline = .true.
  fullestChannelOnly = .true.

  do while (activePipeline)

     localLoop = .true.

     do while (localLoop)
 
        moreItems = .true.

        do while (moreItems)
           call sim_transportItems     ()
           call Pipeline_localProgress ()
           call Pipeline_localGetItems (sim_items (1:sim_itemSize, sim_numItems+1:sim_maxItems), numGetItems)
           sim_numItems = sim_numItems + numGetItems
           moreItems = (sim_numItems > 0)
        end do

        call Pipeline_localFlush    (fullestChannelOnly)
        call Pipeline_localProgress ()
        call Pipeline_localGetItems (sim_items (1:sim_itemSize, sim_numItems+1:sim_maxItems), numGetItems)
        sim_numItems = sim_numItems + numGetItems
        moreItems = (sim_numItems > 0)
        localLoop = moreItems

     end do

     call Pipeline_globalCheckStatus (empty)

     activePipeline = .not.empty

  end do
!
!
!     ...Deactivate/destroy the pipeline.
!
!
  call Pipeline_localDeactivate (doAsyncReturn = .false.)
  call Pipeline_localDestroy    ()
!
!
!     ...Ready!
!
!
  return
end subroutine sim_doPipelineTest
