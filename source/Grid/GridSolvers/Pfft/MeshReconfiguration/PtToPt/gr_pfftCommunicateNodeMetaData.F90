!!****if* source/Grid/GridSolvers/Pfft/MeshReconfiguration/PtToPt/gr_pfftCommunicateNodeMetaData
!!
!! NAME 
!!
!! gr_pfftCommunicateNodeMetaData
!!
!! SYNOPSIS
!!
!! gr_pfftCommunicateNodeMetaData()
!!
!! DESCRIPTION 
!!
!! Communicates the metadata between nodes of FLASH list and nodes of 
!! PFFT list.  The metadata includes send / recv proc IDs so we know which 
!! processor pairs need to communicate.
!!
!! ARGUMENTS
!!
!! SIDE EFFECTS 
!!
!! Adding items to the end of the PFFT list and updating the 
!! attributes of nodes in this list.  We also keep track of the 
!! number of messages we send to each processor in a module level array.
!!
!! NOTES
!! 
!! We must ensure we synchronise all processors before & after we
!! exchange the metadata.  This is because we make use of MPI_Iprobe()
!! during the metadata exchange.  Without synchronisation there is a race
!! condition, and the possibility of intercepting a message originating
!! from another unit. e.g. When we combine PFFT with Multigrid, the order
!! of calls in gr_hgSolve() is gr_hgPfftInitGrid() then
!! gr_hgInitSource().  If a processor reaches gr_hgRestrict() (called by
!! gr_hgInitSource) it will send non-blocking messages which may be
!! intercepted by a processor still in PFFT.  (Tags help us minimise the
!! chance of interception, but the only guarantee is by using barriers)
!!
!!***

subroutine gr_pfftCommunicateNodeMetaData()
#include "constants.h"
#include "Pfft.h"
  use Driver_interface, ONLY : Driver_checkMPIErrorCode, Driver_abortFlash
  use gr_pfftNodeObject, ONLY : node
  use gr_pfftListObject, ONLY : apply_fn_to_nodes
  use gr_pfftData, ONLY : pfft_comm, pfft_myPE, pfft_inLen
  use gr_pfftReconfigData, ONLY : pfft_numFlashNodes, &
       pfft_pencilSize, pfft_logUnit, pfft_listFG, pfft_metaType
  use gr_pfftNodeFnPrototypes, ONLY : gr_pfftFnArgMsgDelivered, &
       gr_pfftCreateNode
  implicit none
#ifdef MANUALLY_PACK_MPI_MESSAGE
  interface
     subroutine PackData(item)
       use gr_pfftNodeObject, ONLY : node
       implicit none
       type(node), pointer :: item
     end subroutine PackData
  end interface
  interface
     subroutine UnPackData(item)
       use gr_pfftNodeObject, ONLY : node
       implicit none
       type(node), pointer :: item
     end subroutine UnPackData
  end interface
#endif
  include "Flash_mpi.h"
  type(node), pointer :: itemRecv, itemSend
  integer, dimension(MPI_STATUS_SIZE) :: statusScalar
  integer, parameter :: pfftTag = 400
  integer :: totPencilSize, recvPencilSize
  integer :: totMsgToSend, fragmentSize
  integer :: numMsgSent, numMsgRecv, numMsgRegistered, ierr
  integer :: srcProc, srcTag, srcMsgSize, srcMsgSizeCheck, destProc
  logical :: msgNotifier, recvMore, sendMore, communicateMore

  if (pfft_myPE == 0) then
#ifdef MANUALLY_PACK_MPI_MESSAGE
     print *, "Manually packing the MPI messages"
#else
     print *, "Using derived data types for the MPI messages"
#endif
  end if


  call MPI_Barrier(pfft_comm(IAXIS), ierr) !*** Important: see warning above subroutine.
  call Driver_checkMPIErrorCode(ierr)


  totMsgToSend = pfft_numFlashNodes
  totPencilSize = pfft_pencilSize  !Lets us know when to break from recv loop.
  recvPencilSize = 0
  numMsgSent = 0 !Number of messages sent (not necessarily delivered).
  numMsgRecv = 0 !Number of messages received (not necessarily delivered).
  numMsgRegistered = 0 !Total of numMsgSent and numMsgRecv.

  sendMore = (totMsgToSend > numMsgSent)
  recvMore = (totPencilSize > recvPencilSize)
  communicateMore = (sendMore .or. recvMore) 

#ifdef DEBUG_PFFT
  write(pfft_logUnit,*) "TotMsgToSend:", totMsgToSend, &
       "TotPencilSize:", totPencilSize
#endif

  do while (communicateMore .eqv. .true.)

     !Let us handle the receives first.
     !--------------------------------------------------------------
     if (recvMore .eqv. .true.) then

        call MPI_Iprobe(MPI_ANY_SOURCE, pfftTag, pfft_comm(IAXIS), &
             msgNotifier, statusScalar, ierr)
        call Driver_checkMPIErrorCode(ierr)


        do while (msgNotifier .eqv. .true.)

           !Extract information about the message matching MPI_Iprobe call.
           srcProc = statusScalar(MPI_SOURCE); srcTag = statusScalar(MPI_TAG)

#ifdef MANUALLY_PACK_MPI_MESSAGE
           call MPI_Get_count(statusScalar, FLASH_INTEGER, srcMsgSize, ierr)
           call Driver_checkMPIErrorCode(pfft_myPE, ierr)
           if (srcMsgSize /= (4+(4*MDIM))) then
              print *, "Processor:", pfft_myPE, "message size:", srcMsgSize
              call Driver_abortFlash &
                   ("[gr_pfftManageNodes]: Unexpected message Size (1)")
           end if 
#else
           call MPI_Get_count(statusScalar, pfft_metaType, srcMsgSize, ierr)
           call Driver_checkMPIErrorCode(ierr)
           if (srcMsgSize /= 1) then
              print *, "Processor:", pfft_myPE, "message size:", srcMsgSize
              call Driver_abortFlash &
                   ("[gr_pfftManageNodes]: Unexpected message Size (1)")
           end if
#endif


#ifdef DEBUG_PFFT
           write(pfft_logUnit,*) "Receiving metadata message from:", srcProc
#endif

           !Allocate space, and then receive the **exact** message 
           !matched by the MPI_IProbe. The itemRecv pointer will always 
           !point to the freshly allocated space.
           call gr_pfftCreateNode(PFFT_PENCIL_NODE, itemRecv)

           !Alternate strategies for receiving the data.
#ifdef MANUALLY_PACK_MPI_MESSAGE
           call MPI_Irecv(itemRecv % packedData(1), (4+(4*MDIM)), FLASH_INTEGER, &
                srcProc, srcTag, pfft_comm(IAXIS), itemRecv % request, ierr)
#else
           call MPI_Irecv(itemRecv % metadata, 1, pfft_metaType, &
                srcProc, srcTag, pfft_comm(IAXIS), itemRecv % request, ierr)
#endif
           call Driver_checkMPIErrorCode(ierr)


           !I need to wait for the message to be delivered now.  The metadata
           !message includes the size of the actual data message that we will
           !later receive.  We accumulate the data message sizes into the  
           !"recvPencilSize" variable.  When this value reaches the size of 
           !the pencil grid assigned to myPE we know that we have received 
           !all metadata messages, and thus we can set "recvMore" to .false..
           call MPI_Wait(itemRecv % request, itemRecv % status, ierr)
           call Driver_checkMPIErrorCode(ierr)
#ifdef DEBUG_PFFT
           write(pfft_logUnit,*) "Metadata message now delivered from:", srcProc
#endif


#ifdef MANUALLY_PACK_MPI_MESSAGE
           call MPI_Get_count(itemRecv % status, FLASH_INTEGER, srcMsgSizeCheck, ierr)
           call Driver_checkMPIErrorCode(ierr)
           if (srcMsgSizeCheck /= (4+(4*MDIM))) then
              print *, "Processor:", pfft_myPE, "message size:", srcMsgSizeCheck
              call Driver_abortFlash &
                   ("[gr_pfftManageNodes]: Unexpected message Size (2)")
           end if
           call UnPackData(itemRecv)  !Place data items in metadata positions.
#else
           call MPI_Get_count(itemRecv % status, pfft_metaType, srcMsgSizeCheck, ierr)
           call Driver_checkMPIErrorCode(ierr)
           if (srcMsgSizeCheck /= 1) then
              print *, "Processor:", pfft_myPE, "message size:", srcMsgSizeCheck
              call Driver_abortFlash &
                   ("[gr_pfftManageNodes]: Unexpected message Size (2)")
           end if
#endif



           !Check that the received pencil coordinates are sensible.
           !In a 2D simulation: pfftStartPos(KAXIS)=pfftEndPos(KAXIS)=1.
           if (any( &
                (itemRecv % metadata % pfftStartPos(1:MDIM) < 1) .or. &
                (itemRecv % metadata % pfftStartPos(1:MDIM) > &
                pfft_inLen(1:MDIM)) .or. &
                (itemRecv % metadata % pfftEndPos(1:MDIM) < 1) .or. &
                (itemRecv % metadata % pfftEndPos(1:MDIM) > &
                pfft_inLen(1:MDIM)) .or. &
                (itemRecv % metadata % pfftEndPos(1:MDIM) < &
                itemRecv % metadata % pfftStartPos(1:MDIM)))) then

              print *, "Processor:", pfft_myPE, &
                   "StartPos:", itemRecv % metadata % pfftStartPos, &
                   "EndPos:", itemRecv % metadata % pfftEndPos, &
                   "Pencil size:", pfft_inLen
              call Driver_abortFlash("Error with received pencil coordinates")
           end if

           fragmentSize = product( &
                (itemRecv % metadata % pfftEndPos(1:MDIM) - &
                itemRecv % metadata % pfftStartPos(1:MDIM) + 1) )
           allocate (itemRecv % buf (fragmentSize))


           !Now we have finished with our pointer we set it to null.
           !We do not deallocate the memory because the received message 
           !will sit in this space.  We can retrieve this data later by 
           !following the head and tail list pointers.
           nullify(itemRecv)


           !Advance the message counters.
           numMsgRecv = numMsgRecv + 1
           numMsgRegistered = numMsgRegistered + 1
           recvPencilSize = recvPencilSize + fragmentSize


           !Test if there is another message available.
           call MPI_Iprobe(MPI_ANY_SOURCE, pfftTag, pfft_comm(IAXIS), &
                msgNotifier, statusScalar, ierr)
           call Driver_checkMPIErrorCode(ierr)
        end do

        recvMore = (totPencilSize > recvPencilSize)

     end if



     !Let us now handle one send.  After the send we will attempt to receive
     !more messages.  Iterate through the send buffer and send an individual 
     !sections to the relevant process.
     !--------------------------------------------------------------
     if (sendMore .eqv. .true.) then

        if (numMsgSent == 0) then
           itemSend => pfft_listFG % H  !Point to start of FLASH list.
        end if

        !Check that this node is usable:
        if (.not.associated(itemSend)) then
           call Driver_abortFlash("Item not associated!!!")
        end if

        destProc = itemSend % metadata % pfftProcID

#ifdef DEBUG_PFFT
        write(pfft_logUnit,*) "Sending metadata message to proc:", destProc
#endif

        !Alternate strategies for sending the data.
#ifdef MANUALLY_PACK_MPI_MESSAGE
        call PackData(itemSend)
        call MPI_Isend(itemSend % packedData(1), (4+(4*MDIM)), FLASH_INTEGER, &
             destProc, pfftTag, pfft_comm(IAXIS), itemSend % request, ierr)
#else
        call MPI_Isend(itemSend % metadata, 1, pfft_metaType, &
             destProc, PfftTag, pfft_comm(IAXIS), itemSend % request, ierr)
#endif
        call Driver_checkMPIErrorCode(ierr)


        !Increment message counters.
        numMsgSent = numMsgSent + 1
        numMsgRegistered = numMsgRegistered + 1
        itemSend => itemSend % next

        sendMore = (numMsgSent < totMsgToSend)
     end if   !Testing for messages to send.
     !-------------------------------------------------------------- 

     communicateMore = (sendMore .or. recvMore) 

  end do   !Loop until we account for all messages involving myPE.


  if (totPencilSize /= recvPencilSize) then
     print *, "Mismatch: totPencilSize=", totPencilSize, &
          "recvPencilSize=", recvPencilSize
     call Driver_abortFlash("Too much metadata on one processor")
  end if

  !We only need to wait for delivery of the messages from the FLASH list.
  !The earlier MPI_Wait in the "recvMore" loop has already guaranteed 
  !delievery of PFFT list messages.
  call apply_fn_to_nodes(gr_pfftFnArgMsgDelivered, pfft_listFG)


  call MPI_Barrier(pfft_comm(IAXIS), ierr) !*** Important: see warning above subroutine.
  call Driver_checkMPIErrorCode(ierr)

end subroutine gr_pfftCommunicateNodeMetaData


#ifdef MANUALLY_PACK_MPI_MESSAGE
subroutine PackData(item)
  use gr_pfftNodeObject, ONLY : node
  implicit none
  type(node), pointer :: item
  item % packedData(1) = item % metadata % flashProcID
  item % packedData(2) = item % metadata % flashBlockID
  item % packedData(3) = item % metadata % pfftProcID
  item % packedData(4) = item % metadata % tagID
  item % packedData(5:7) = item % metadata % flashStartPos(1:MDIM)
  item % packedData(8:10) = item % metadata % flashEndPos(1:MDIM)
  item % packedData(11:13) = item % metadata % pfftStartPos(1:MDIM)
  item % packedData(14:16) = item % metadata % pfftEndPos(1:MDIM)
end subroutine PackData


subroutine UnPackData(item)
  use gr_pfftNodeObject, ONLY : node
  implicit none
  type(node), pointer :: item
  item % metadata % flashProcID = item % packedData(1)
  item % metadata % flashBlockID = item % packedData(2)
  item % metadata % pfftProcID = item % packedData(3)
  item % metadata % tagID = item % packedData(4)
  item % metadata % flashStartPos(1:MDIM) = item % packedData(5:7)
  item % metadata % flashEndPos(1:MDIM) = item % packedData(8:10)
  item % metadata % pfftStartPos(1:MDIM) = item % packedData(11:13)
  item % metadata % pfftEndPos(1:MDIM) = item % packedData(14:16)
end subroutine UnPackData
#endif
