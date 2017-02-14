!!****f* source/flashUtilities/Pipeline/Pipeline_localSendItem
!!
!! NAME
!!  
!!  Pipeline_localSendItem
!!
!! SYNOPSIS
!! 
!!  call Pipeline_localSendItem (real,    intent (in)  :: userItem (:),
!!                               integer, intent (in)  :: userProcID,
!!                               logical, intent (out) :: isHandled)
!!
!! DESCRIPTION
!!
!!  The routine will try to add the user supplied item to the local send buffer on
!!  the pipeline channel corresponding to the user supplied channel processor procID.
!!  This can be done, if there is no pending send on the send buffer and if there
!!  is enough space on the send buffer. This part of the routine can be considered
!!  the local feeder of the pipeline.
!!
!! ARGUMENTS
!!
!!  userItem   : the item (array of elements) to be added to the send buffer
!!  userProcID : channel processor ID
!!  isHandled  : is true, if the item was successfully added to the buffer
!!
!! NOTES
!!
!!  If the item supplied by the user contains more elements than the pipeline items
!!  can handle, the program aborts.
!!
!!  It is the user's responsibility to make sure that the supplied channel processor
!!  ID actually belongs to a channel of the pipeline on the local processor. If the
!!  userProcID channel is not found within all the channels of the pipeline on the
!!  local processor, then currently we signal that something must have gone wrong.
!!  Either the pipeline for the local processor was set up incorrectly or userProcID
!!  does not correspond to a channel processor. In this case the only solution is to
!!  abort the program.
!!
!!  In the future we might just simply ignore the sending item request, if userProcID
!!  is not one of the channels. We simply return isHandled = false and continue with
!!  the program. 
!!
!!  Currently the routine is written in such a way that the item will always be added
!!  to the send buffer, i.e. isHandled should always be true. In case there is a
!!  pending send of the buffer, we call the progress sending communication repeatedly
!!  until the send for the particular channel has completed (i.e. the send buffer is
!!  free to be reused). A return of isHandled equal to false indicates an abnormality
!!  and should be catched by the calling routine.
!!
!!***

subroutine Pipeline_localSendItem (userItem, userProcID,   isHandled)

  implicit none

  real,    intent (in)  :: userItem (:)
  integer, intent (in)  :: userProcID
  logical, intent (out) :: isHandled

  isHandled = .false.

  return
end subroutine Pipeline_localSendItem
