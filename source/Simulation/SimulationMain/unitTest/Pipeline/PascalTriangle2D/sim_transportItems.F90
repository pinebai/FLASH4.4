!!****if* source/Simulation/SimulationMain/unitTest/Pipeline/PascalTriangle2D/sim_transportItems
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
!!   This routine transports the items through the pipeline on the local processor.
!!   The rule of transportation depends on the local processor rank and is such
!!   that it is down the pipeline tree:
!!
!!                                         0
!!                                        / \
!!                                       1   2
!!                                      / \ / \
!!                                     3   4   5
!!                                    / \ / \ / \
!!                                   6   7   8   9
!!                                  / \ / \ / \ / \
!!                                10  11  12  13  14
!!
!!   Note that on a given rank there are either 2 sending channels (in the example
!!   picture ranks 0 to 9) or no sending channels (in the example picture ranks 10
!!   to 14). When 2 sending channels are present and there are n items, we send the
!!   first n/2 items through the lowest rank channel and the last n/2 items through
!!   the highest rank channel. If n is odd, we will send int(n/2)+1 items through
!!   the lowest rank channel and int(n/2) items through the highest rank channel.
!!   When no sending channels are present, the items are 'removed' by resetting the
!!   item count to zero.
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
  integer :: lowRank, highRank
  integer :: nLow
  integer :: numSend
!
!
!     ...Select the rank and act accordingly.
!
!
  nLow = int (sim_numItems / 2) + mod (sim_numItems,2)
  numSend = 0

  if (sim_numSendingChannels == 2) then

      lowRank = sim_sendingChannels (1)

      do item = 1, nLow
         call Pipeline_localSendItem (sim_items (1:sim_itemSize,item), lowRank, isHandled)
         if (.not.isHandled) then
              call Driver_abortFlash ('[sim_transportItems] ERROR: Item send failure!')
         end if
         numSend = numSend + 1
      end do

      highRank = sim_sendingChannels (2)

      do item = nLow+1, sim_numItems
         call Pipeline_localSendItem (sim_items (1:sim_itemSize,item), highRank, isHandled)
         if (.not.isHandled) then
              call Driver_abortFlash ('[sim_transportItems] ERROR: Item send failure!')
         end if
         numSend = numSend + 1
      end do

  else if (sim_numSendingChannels == 0) then

      numSend = sim_numItems
      sim_storedItems = sim_storedItems + numSend

  else
      call Driver_abortFlash ('[sim_transportItems] ERROR: # of sending channels .ne. 0 or 2!')
  end if

  sim_numItems = sim_numItems - numSend
!
!
!     ...Ready!
!
!
  return
end subroutine sim_transportItems
