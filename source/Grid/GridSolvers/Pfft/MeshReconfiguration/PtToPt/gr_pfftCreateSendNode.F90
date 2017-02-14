!!****if* source/Grid/GridSolvers/Pfft/MeshReconfiguration/PtToPt/gr_pfftCreateSendNode
!!
!! NAME 
!!
!! gr_pfftCreateSendNode
!!
!! SYNOPSIS
!!
!! gr_pfftCreateSendNode(integer (IN) :: flashProcID, &
!!                       integer (IN) :: flashBlockID, &
!!                       integer (IN) :: flashStartPos(MDIM), &
!!                       integer (IN) :: flashEndPos(MDIM), &
!!                       integer (IN) :: pfftProcID, &
!!                       integer (IN) :: pfftStartPos(MDIM), &
!!                       integer (IN) :: pfftEndPos(MDIM))
!!
!! DESCRIPTION 
!!
!! Creates a new send node which we place on FLASH list by calling 
!! gr_pfftCreateNode.
!!
!! ARGUMENTS
!!
!! flashProcID - My processor ID
!! flashBlockID - Block ID for this fragment.
!! flashStartPos - X,Y,Z start index of fragment within FLASH block.
!! flashEndPos - X,Y,Z end index of fragment within FLASH block.
!! pfftProcID - PFFT processor due to receive this fragment.
!! pfftStartPos - X,Y,Z start index of fragment within PFFT block.
!! pfftEndPos - X,Y,Z end index of fragment within PFFT block.
!!
!! SIDE EFFECTS 
!!
!! Adding items to the end of the FLASH list and updating the 
!! attributes of nodes in this list.  We also keep track of the 
!! number of messages we send to each processor in a module level array.
!!
!! NOTES
!! 
!! We maintain two lists: one list is used to hold FLASH block fragments, 
!! and the other list is used to hold PFFT block fragments.  The lists are 
!! retained in memory so that the same communication pattern (stored in the
!! lists) can be re-used for multiple parallel FFTs.
!!
!! We can send more than 1 fragment to a single processor because we
!! maintain a tagID within each node.  For simplicity we just perform
!! an MPI communication even if sendPE = recvPE.
!!
!!***

subroutine gr_pfftCreateSendNode(flashProcID, flashBlockID, flashStartPos, &
     flashEndPos, pfftProcID, pfftStartPos, pfftEndPos)
#include "constants.h"
#include "Pfft.h"
  use Driver_interface, ONLY : Driver_abortFlash
  use gr_pfftNodeObject, ONLY : node
  use gr_pfftNodeFnPrototypes, ONLY : gr_pfftCreateNode
  use gr_pfftReconfigData, ONLY : pfft_numMsgSendToEachProc
  implicit none
  integer, intent(IN) :: flashProcID, flashBlockID
  integer, dimension(1:MDIM), intent(IN) :: flashStartPos, flashEndPos
  integer, intent(IN) :: pfftProcID
  integer, dimension(1:MDIM), intent(IN) :: pfftStartPos, pfftEndPos
  type(node), pointer :: item
  integer :: flashFragmentSize, pfftFragmentSize

  flashFragmentSize = &
       product( (flashEndPos(1:MDIM) - flashStartPos(1:MDIM) + 1) )

  pfftFragmentSize = &
       product( (pfftEndPos(1:MDIM) - pfftStartPos(1:MDIM) + 1) )

  if (flashFragmentSize /= pfftFragmentSize) then
     call Driver_abortFlash("Mismatch in fragment calculation!")
  end if

  !Create a node and add it to the FLASH list.
  call gr_pfftCreateNode(PFFT_FLASH_NODE, item)

  item % metadata % flashProcID = flashProcID
  item % metadata % flashBlockID = flashBlockID
  item % metadata % flashStartPos = flashStartPos
  item % metadata % flashEndPos = flashEndPos
  item % metadata % pfftProcID = pfftProcID
  item % metadata % pfftStartPos = pfftStartPos
  item % metadata % pfftEndPos = pfftEndPos

  pfft_numMsgSendToEachProc(pfftProcID) = &
       pfft_numMsgSendToEachProc(pfftProcID) + 1

  item % metadata % tagID = pfft_numMsgSendToEachProc(pfftProcID)
  allocate(item % buf (flashFragmentSize))

end subroutine gr_pfftCreateSendNode
