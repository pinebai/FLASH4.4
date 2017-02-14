!!****if* source/Simulation/SimulationMain/unitTest/Pipeline/Linear/sim_transportItems
!!
!!  NAME 
!!
!!   sim_transportItems
!!
!!  SYNOPSIS
!!
!!   sim_transportItems ()
!!
!!  DESCRIPTION
!!
!!   This routine transports the items through the linear pipeline on the local
!!   processor. In a linear pipeline each processor has exactly one sending channel
!!   except for the processor with the highest rank, which is only responsible for
!!   removing (resetting the item count to zero) the items from the pipeline. 
!!
!! ARGUMENTS
!!
!!***

subroutine sim_transportItems ()

  use Driver_interface,    ONLY : Driver_abortFlash

  use Pipeline_interface,  ONLY : Pipeline_localSendItem

  use Simulation_data,     ONLY : sim_items,              &
                                  sim_itemSize,           &
                                  sim_numItems,           &
                                  sim_numSendingChannels, &
                                  sim_sendingChannels,    &
                                  sim_storedItems

  implicit none

  logical :: isHandled

  integer :: item
  integer :: numSend
  integer :: sendRank
!
!
!     ...Select the rank and act accordingly.
!
!
  numSend = 0

  if (sim_numSendingChannels == 1) then

      sendRank = sim_sendingChannels (1)

      do item = 1, sim_numItems
         call Pipeline_localSendItem (sim_items (1:sim_itemSize,item), sendRank, isHandled)
         if (.not.isHandled) then
              call Driver_abortFlash ('[sim_transportItems] ERROR: Item send failure!')
         end if
         numSend = numSend + 1
      end do

  else if (sim_numSendingChannels == 0) then

      numSend = sim_numItems
      sim_storedItems = sim_storedItems + numSend

  else
      call Driver_abortFlash ('[sim_transportItems] ERROR: # of sending channels .ne. 0 or 1!')
  end if

  sim_numItems = sim_numItems - numSend
!
!
!     ...Ready!
!
!
  return
end subroutine sim_transportItems
