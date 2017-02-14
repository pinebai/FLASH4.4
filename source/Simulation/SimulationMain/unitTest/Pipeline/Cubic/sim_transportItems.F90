!!****if* source/Simulation/SimulationMain/unitTest/Pipeline/Cubic/sim_transportItems
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
!!   that it is through the cube. On a given rank there can be up to 3 sending channels.
!!   The following action is taken for n items:
!!
!!   no sending channel: all n items are 'removed' by resetting the item
!!                       count to zero
!!
!!   1 sending channel : send all n items throught that channel
!!
!!   2 sending channels: lowest  channel -> send int[(n+1)/2] items
!!                       highest channel -> send int[ n   /2] items
!!
!!   3 sending channels: lowest  channel -> send int[(n+2)/3] items
!!                       middle  channel -> send int[(n+1)/3] items
!!                       highest channel -> send int[ n   /3] items
!!
!!   The sending channels have already been determined in lowest, middle, highest order.
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
  integer :: nLow, nMid, nHgh
  integer :: numSend
  integer :: rank, lowRank, midRank, hghRank
!
!
!     ...Select the cases.
!
!
  numSend = 0

  select case (sim_numSendingChannels)

  case (3)

    nLow = (sim_numItems + 2) / 3
    nMid = (sim_numItems + 1) / 3
    nHgh =  sim_numItems      / 3

    if (nLow + nMid + nHgh /= sim_numItems) then
        call Driver_abortFlash ('[sim_transportItems] ERROR: Mismatch # of sending items!')
    end if

    lowRank = sim_sendingChannels (1)
    midRank = sim_sendingChannels (2)
    hghRank = sim_sendingChannels (3)

    do item = 1, nLow
       call Pipeline_localSendItem (sim_items (1:sim_itemSize,item), lowRank, isHandled)
       if (.not.isHandled) then
            call Driver_abortFlash ('[sim_transportItems] ERROR: Item send failure!')
       end if
       numSend = numSend + 1
    end do

    do item = nLow + 1, nLow + nMid
       call Pipeline_localSendItem (sim_items (1:sim_itemSize,item), midRank, isHandled)
       if (.not.isHandled) then
            call Driver_abortFlash ('[sim_transportItems] ERROR: Item send failure!')
       end if
       numSend = numSend + 1
    end do

    do item = nLow + nMid + 1, sim_numItems
       call Pipeline_localSendItem (sim_items (1:sim_itemSize,item), hghRank, isHandled)
       if (.not.isHandled) then
            call Driver_abortFlash ('[sim_transportItems] ERROR: Item send failure!')
       end if
       numSend = numSend + 1
    end do

  case (2)

    nLow = (sim_numItems + 1) / 2
    nHgh =  sim_numItems      / 2

    if (nLow + nHgh /= sim_numItems) then
        call Driver_abortFlash ('[sim_transportItems] ERROR: Mismatch # of sending items!')
    end if

    lowRank = sim_sendingChannels (1)
    hghRank = sim_sendingChannels (2)

    do item = 1, nLow
       call Pipeline_localSendItem (sim_items (1:sim_itemSize,item), lowRank, isHandled)
       if (.not.isHandled) then
            call Driver_abortFlash ('[sim_transportItems] ERROR: Item send failure!')
       end if
       numSend = numSend + 1
    end do

    do item = nLow + 1, sim_numItems
       call Pipeline_localSendItem (sim_items (1:sim_itemSize,item), hghRank, isHandled)
       if (.not.isHandled) then
            call Driver_abortFlash ('[sim_transportItems] ERROR: Item send failure!')
       end if
       numSend = numSend + 1
    end do

  case (1)

    rank = sim_sendingChannels (1)

    do item = 1, sim_numItems
       call Pipeline_localSendItem (sim_items (1:sim_itemSize,item), rank, isHandled)
       if (.not.isHandled) then
            call Driver_abortFlash ('[sim_transportItems] ERROR: Item send failure!')
       end if
       numSend = numSend + 1
    end do

  case (0)

    numSend = sim_numItems
    sim_storedItems = sim_storedItems + numSend

  case default
       call Driver_abortFlash ('[sim_transportItems] ERROR: # of sending channels .ne. 0,1,2,3!')
  end select

  sim_numItems = sim_numItems - numSend
!
!
!     ...Ready!
!
!
  return
end subroutine sim_transportItems
