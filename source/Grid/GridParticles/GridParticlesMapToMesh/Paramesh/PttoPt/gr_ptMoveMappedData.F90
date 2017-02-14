!!****if* source/Grid/GridParticles/GridParticlesMapToMesh/Paramesh/PttoPt/gr_ptMoveMappedData
!!
!! NAME
!!  gr_ptMoveMappedData
!!
!! SYNOPSIS
!!
!!  gr_ptMoveMappedData(integer,intent(IN) :: varGrid, &
!!                       integer,intent(IN) :: bufferSize, &
!!                       real,dimension(bufferSize),intent(INOUT) :: sendBuf, &
!!                       integer,intent(INOUT) :: sendCount, &
!!                       real,dimension(bufferSize),intent(INOUT) :: recvBuf)
!!
!! DESCRIPTION
!!
!! Routine which manages the communication of smeared grid cells between processes. 
!! The smeared grid cells are stored in different sections of a send buffer, and each section 
!! is sent to the appropriate process.  This implementation uses non-blocking point to point 
!! messaging with no global communication.  We refer to "gr_ptNumMessagesToSend" for the 
!! number of messages we need to send, and "gr_ptRecvSpecifier" for the number of messages 
!! we need to receive.  We do not know who we will receive from so we use 
!! MPI_ANY_SOURCE to pick up any messages that have been sent. 
!! 
!!
!! ARGUMENTS
!!               varGrid:   Index of gridded variable to receive interpolated
!!                              quantity
!!               bufferSize:  The size of the sendBuf and recvBuf arrays
!!               sendBuf:  An array used to store data intended for another processor
!!               sendCount:  The number of data elementes to be sent to another
!!                           processor
!!               recvBuf:  An array containing the data just receieved from another 
!!                         processor
!!
!! PARAMETERS
!! 
!!***

subroutine gr_ptMoveMappedData(varGrid,bufferSize,sendBuf,sendCount,recvBuf)

  use Grid_data, ONLY : gr_meshMe, gr_meshNumProcs, gr_meshComm
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Logfile_interface, ONLY: Logfile_stampMessage
  use Driver_interface, ONLY : Driver_abortFlash
  use gr_ptInterface, ONLY : gr_ptPackUnpackData  
  use gr_ptMapData, ONLY : gr_ptRecvSpecifier, gr_ptNumMessagesToSend

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
#include "gr_ptMapToMesh.h"

  integer,intent(IN) :: varGrid
  integer,intent(IN) :: bufferSize
  real,dimension(bufferSize),intent(INOUT) :: sendBuf
  integer,intent(INOUT) :: sendCount
  real,dimension(bufferSize),intent(INOUT) :: recvBuf

  integer, parameter :: MapTag = 325
  integer, dimension(LOW:HIGH,MDIM) :: sectionCoords
  integer, dimension(MPI_STATUS_SIZE) :: statusScalar   
  integer, allocatable, dimension(:,:) :: statusArray
  integer, allocatable, dimension(:) :: requestArray

  integer :: totMsgToSend, totMsgToRecv, totMsgInvolvingMyPE
  integer :: numMsgSent, numMsgRecv, numMsgRegistered 
  integer :: recvBufPtr, sendBufPtr, headerPtr, error, ierr
  integer :: srcProc, srcTag, msgID, srcMsgSize, srcMsgSizeCheck
  integer :: destProc, elementsInMessage, numbElements, recvCount
  integer :: lastMsgSentAndDelivered, lastMsgSentID
  logical :: msgNotifier, msgDelivered

  !Switch on both message throttling modes for safety.
#define THROTTLE_TO_ONE_RECV
#define THROTTLE_TO_ONE_SEND


  if (gr_meshNumProcs == 1) then
     return   !No communication necessary.
  else
     !We use MPI_IProbe later, so we must ensure we only probe for messages 
     !from within this subroutine.  This means we must insert a barrier
     !at the start and end of this subroutine.
     call MPI_Barrier(gr_meshComm, ierr)
  end if

  !Determine the total number of messages we are destined to send and receive.
  totMsgToSend = gr_ptNumMessagesToSend
  totMsgToRecv = gr_ptRecvSpecifier(gr_meshMe)
  totMsgInvolvingMyPE = totMsgToSend + totMsgToRecv


  !Return if we do not participate in the forthcoming messaging madness.
  !Everyone must reach the MPI_BARRIER at the end before leaving!
  if (totMsgInvolvingMyPE > 0) then

     !Size status and request arrays appropriately.
     allocate(statusArray(MPI_STATUS_SIZE, totMsgInvolvingMyPE), &
          requestArray(totMsgInvolvingMyPE), STAT=error)
     if (error /= 0) then
        call Driver_abortFlash("[gr_ptMoveMappedData]: status/request arrays cannot be allocated!")
     end if

     !Necessary initialisations:
     numMsgSent = 0 !Number of messages we have sent (not necessarily delivered).
     numMsgRecv = 0 !Number of messages we have received (not necessarily delivered).
     numMsgRegistered = 0 !Total of numMsgSent and numMsgRecv.
     lastMsgSentID = 0 !Request object ID of the last message we sent.
     lastMsgSentAndDelivered = 0 !Last message we sent which has been delivered.
     recvBufPtr = 1 !Integer counter in the receive buffer.
     sendBufPtr = 1 !Integer counter in the send buffer.
     recvCount = 0 !Total number of FLASH_REALs we have received.


     !Remain in send/recv loop until all messages are accounted for.
     !--------------------------------------------------------------
     do while (numMsgRegistered < totMsgInvolvingMyPE)

115     continue  !! This is the start of the communication loop.


        !Let us handle the receives first.
        !--------------------------------------------------------------
        if (numMsgRecv < totMsgToRecv) then

           call MPI_Iprobe(MPI_ANY_SOURCE, MapTag, gr_meshComm, &
                msgNotifier, statusScalar, ierr)

           do while (msgNotifier .eqv. .true.)

              !Extract information about the message matching the MPI_Iprobe pattern.
              srcProc = statusScalar(MPI_SOURCE); srcTag = statusScalar(MPI_TAG)
              call MPI_Get_count(statusScalar, FLASH_REAL, srcMsgSize, ierr)


              !Receive the **exact** message matched by the MPI_IProbe.
              !This is essential because the point to point messages are not all the same size. 
              msgID = numMsgRegistered + 1
              call MPI_Irecv(recvBuf(recvBufPtr), srcMsgSize, FLASH_REAL, &
                   srcProc, srcTag, gr_meshComm, requestArray(msgID), ierr)


#ifdef THROTTLE_TO_ONE_RECV
              !This should be switched on!  I'm unsure if the next MPI_Iprobe can probe the 
              !same message if the message has yet to be delivered into myPE's recvBuf.  Issuing 
              !an extra MPI_Wait() ensures the message is delivered before the next MPI_Iprobe!
              call MPI_Wait(requestArray(msgID), statusScalar, ierr)
              call MPI_Get_count(statusScalar, FLASH_REAL, srcMsgSizeCheck, ierr)
              if (srcMsgSizeCheck /= srcMsgSize) then
                 print *, "We have somehow received a 2nd message from processor:", &
                      srcProc, "before receiving the 1st message from the same processor"
                 call Driver_abortFlash & 
                      ("[gr_ptMoveMappedData]: Message ordering problem - MPI issue!")
              endif
#endif THROTTLE_TO_ONE_RECV


              !Advance the receive buffer pointer and the message counters.
              recvBufPtr = recvBufPtr + srcMsgSize
              numMsgRecv = numMsgRecv + 1
              numMsgRegistered = numMsgRegistered + 1
              recvCount = recvCount + srcMsgSize  !Required for gr_ptPackUnpackData().


              !Test if there is another message available.
              call MPI_Iprobe(MPI_ANY_SOURCE, MapTag, gr_meshComm, &
                   msgNotifier, statusScalar, ierr)

           end do
        end if   !Testing for messages to receive.
        !--------------------------------------------------------------           


        !Let us now handle one send.  After the send we will attempt to receive
        !more messages.  Iterate through the send buffer and send an individual 
        !sections to the relevant process.
        !--------------------------------------------------------------
        if (numMsgSent < totMsgToSend) then

#ifdef THROTTLE_TO_ONE_SEND
           !May be required in case we ever issue so many sends that we 
           !overflow the sending process' system buffer.
           if (lastMsgSentAndDelivered < numMsgSent) then
              !Test the request object and if the message is not yet delivered 
              !go back to the main communication loop.  While we wait for the condition to be 
              !true: attempt to receive more messages, or if there are no more messages 
              !to receive then we just busy wait.  Using an MPI_Wait in this position in the code 
              !could lead to deadlock, so we use the non-blocking MPI_Test.
              call MPI_Test(requestArray(lastMsgSentID), msgDelivered, statusScalar, ierr)
              if (msgDelivered .eqv. .false.) then
                 goto 115   !! Warning.  Go to start of communication loop
              else
                 lastMsgSentAndDelivered = lastMsgSentAndDelivered + 1
                 if (lastMsgSentAndDelivered /= numMsgSent) & 
                      call Driver_abortFlash("[gr_ptMoveMappedData]: Problem with message counters!")
              end if
           end if
#endif THROTTLE_TO_ONE_SEND


           headerPtr = sendBufPtr

           !Each message consists of metadata and data.
           destProc = int(sendBuf(headerPtr+BLKPROC-1))
           sectionCoords(LOW,IAXIS) = int(sendBuf(headerPtr+COORDSID-1))
           sectionCoords(HIGH,IAXIS) = int(sendBuf(headerPtr+COORDSID))
           sectionCoords(LOW,JAXIS) = int(sendBuf(headerPtr+COORDSID+1))
           sectionCoords(HIGH,JAXIS) = int(sendBuf(headerPtr+COORDSID+2))
           sectionCoords(LOW,KAXIS) = int(sendBuf(headerPtr+COORDSID+3))
           sectionCoords(HIGH,KAXIS) = int(sendBuf(headerPtr+COORDSID+4))

           numbElements = (sectionCoords(HIGH,IAXIS)-sectionCoords(LOW,IAXIS)+1) * &
                (sectionCoords(HIGH,JAXIS)-sectionCoords(LOW,JAXIS)+1) * &
                (sectionCoords(HIGH,KAXIS)-sectionCoords(LOW,KAXIS)+1)

           if (destProc == NONEXISTENT) then
              !This should not happen as gr_ptFindNegh now uses gr_getBlkHandle() 
              !which queries PARAMESH cached data to find all neighbors (block & processor).
              call Driver_abortFlash("[gr_ptMoveMappedData]: Destination processor unknown!")
           end if
           elementsInMessage = SIZE_HEADER + numbElements
           msgID = numMsgRegistered + 1


           !Send message and record the request object ID.  We use this request 
           !object ID to guarentee delivery of this message (required 
           !in THROTTLE_TO_ONE_SEND mode).
           call MPI_Isend(sendBuf(sendBufPtr), elementsInMessage, FLASH_REAL, &
                destProc, MapTag, gr_meshComm, requestArray(msgID), ierr)
           lastMsgSentID = msgID


           !Advance the send buffer pointer and increment message counters.
           sendBufPtr = sendBufPtr + elementsInMessage
           numMsgSent = numMsgSent + 1
           numMsgRegistered = numMsgRegistered + 1


           !Ensure we stay within the bounds of valid data in the send buffer.
           if((sendBufPtr-1) > sendCount) then
              call Driver_abortFlash("[gr_ptPackUnpackData] Overrun send buffer.") 
           end if

        end if   !Testing for messages to send.
        !-------------------------------------------------------------- 

     end do   !Loop until we account for all messages involving myPE.


     !We just need to wait for successful delivery now.
     call MPI_Waitall(totMsgInvolvingMyPE, requestArray, statusArray, ierr)


     !Now we can be assured that all communication is complete.
     !Use the same function from the MoveSieve implementation to unpack data.
     !NOTE: We will not repack any data in this implementation.
     call gr_ptPackUnpackData(varGrid, bufferSize, sendBuf, sendCount, recvBuf, recvCount)


     deallocate(statusArray, requestArray, STAT=error)
     if (error /= 0) then
        call Driver_abortFlash("[gr_ptMoveMappedData]: status/request arrays cannot be deallocated!")
     end if

  end if

  !We use MPI_IProbe previously, so we must ensure we only probe for messages
  !from within this subroutine.  This means we must insert a barrier
  !at the start and end of this subroutine.
  call MPI_Barrier(gr_meshComm, ierr)

end subroutine gr_ptMoveMappedData
