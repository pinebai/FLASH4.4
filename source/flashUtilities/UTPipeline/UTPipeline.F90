#define ASSERT(x) \
  if (.not.x) call Driver_abortFlash("Assertion fail")
#define ASSERT_MPI_SUCCESS(ierr) \
  call Driver_checkMPIErrorCode(ierr)

module UTPipeline

#ifdef UTPIPELINE_UNIT_TEST
  use mpi
#else
  use Driver_interface, ONLY : Driver_abortFlash, Driver_checkMPIErrorCode
#endif

  implicit none

#ifdef UTPIPELINE_UNIT_TEST
  integer, parameter :: FLASH_INTEGER = MPI_INTEGER
  integer, parameter :: FLASH_REAL = MPI_DOUBLE_PRECISION !WARNING MUST PROMOTE REALS TO DP
#else
  include 'Flash_mpi.h'
#endif

  real, allocatable, save, dimension(:,:,:) :: utpipe_sendBuf
  real, allocatable, save, dimension(:,:,:) :: utpipe_recvBuf
  real, allocatable, save, dimension(:,:) :: utpipe_itemBuf

  integer, parameter :: utpipe_tag = 1235
  integer, allocatable, save, dimension(:,:) :: utpipe_recvStatus
  integer, allocatable, save, dimension(:) :: utpipe_recvRequest
  integer, allocatable, save, dimension(:) :: utpipe_recvIndex
  integer, allocatable, save, dimension(:) :: utpipe_recvCount

  integer, allocatable, save, dimension(:,:) :: utpipe_sendStatus
  integer, allocatable, save, dimension(:) :: utpipe_sendRequest
  integer, allocatable, save, dimension(:) :: utpipe_sendIndex
  integer, allocatable, save, dimension(:) :: utpipe_sendCount

  integer, allocatable, save, dimension(:) :: utpipe_procList
  integer, save :: utpipe_itemCount

  integer, save :: utpipe_comm
  integer, save :: utpipe_size
  integer, save :: utpipe_rank

  integer, save :: utpipe_itemSize
  integer, save :: utpipe_maxItems
  integer, save :: utpipe_channelSize
  integer, save :: utpipe_numChannels

  logical, save :: utpipe_isSendCommDone
  logical, save :: utpipe_isRecvCommDone

  logical, save :: utpipe_isInitialized = .false.
  logical, save :: utpipe_isCommInitialized = .false.
  logical, save :: utpipe_isCommDone = .false.

  integer, save :: utpipe_logUnit
  logical, save :: utpipe_doLog

  integer, allocatable, save, dimension(:) :: utpipe_sendState
  integer, parameter :: OPEN_STATE = -3000
  integer, parameter :: PROMISE_TO_CLOSE_STATE = -4000
  integer, parameter :: WAITING_TO_CLOSE_STATE = -5000
  integer, parameter :: CLOSE_STATE = -6000

contains

  !Needs to be called whenever the FLASH grid changes or whenever
  !the pipeline code is used for another purpose.
  subroutine UTPipeline_init(itemSize, maxItems, channelSize, comm, numChannels, procList, logUnit)
    implicit none
    integer, intent(IN) :: itemSize, maxItems, channelSize, comm, numChannels
    integer, dimension(:), intent(IN) :: procList
    integer, optional, intent(IN) :: logUnit
    integer :: ierr

    utpipe_itemSize = itemSize
    utpipe_maxItems = maxItems
    utpipe_channelSize = channelSize
    utpipe_comm = comm
    utpipe_numChannels = numChannels
    
    utpipe_logUnit = -1
    utpipe_doLog = .false.
    if (present(logUnit)) then
       utpipe_doLog = (logUnit >= 0) !Must be a valid Fortran unit
       if (utpipe_doLog) utpipe_logUnit = logUnit
    end if

    if (numChannels > 0) then
       allocate(utpipe_itemBuf(itemSize,maxItems))
       allocate(utpipe_sendBuf(itemSize,channelSize,numChannels))
       allocate(utpipe_sendStatus(MPI_STATUS_SIZE,numChannels))
       allocate(utpipe_sendRequest(numChannels))
       allocate(utpipe_sendIndex(numChannels))
       allocate(utpipe_sendCount(numChannels))
       allocate(utpipe_sendState(numChannels))
       allocate(utpipe_recvBuf(itemSize,channelSize,numChannels))
       allocate(utpipe_recvStatus(MPI_STATUS_SIZE,numChannels))
       allocate(utpipe_recvRequest(numChannels))
       allocate(utpipe_recvIndex(numChannels))
       allocate(utpipe_recvCount(numChannels))
       allocate(utpipe_procList(numChannels))
       utpipe_procList(1:numChannels) = procList(1:numChannels)
       utpipe_sendRequest = MPI_REQUEST_NULL
       utpipe_recvRequest = MPI_REQUEST_NULL
    end if

    call MPI_Comm_rank(utpipe_comm, utpipe_rank, ierr)
    ASSERT_MPI_SUCCESS(ierr)
    call MPI_Comm_size(utpipe_comm, utpipe_size, ierr)
    ASSERT_MPI_SUCCESS(ierr)

    utpipe_isInitialized = .true.
  end subroutine UTPipeline_init


  subroutine UTPipeline_initComm()
    implicit none
    integer :: i
    if (.not.utpipe_isInitialized) then
       call Driver_abortFlash("Must intialize the pipeline first")
    end if

    utpipe_itemCount = 0
    if (utpipe_numChannels > 0) then
       utpipe_sendCount(:) = 0
       utpipe_recvCount(:) = 0
       utpipe_isCommDone = .false.
       utpipe_isSendCommDone = .false.
       utpipe_isRecvCommDone = .false.
       do i = 1, utpipe_numChannels
          utpipe_sendState(i) = OPEN_STATE
          call utpipe_postRecvMsg(i)
       end do
    else
       utpipe_isCommDone = .true.
       utpipe_isSendCommDone = .true.
       utpipe_isRecvCommDone = .true.
    end if
    utpipe_isCommInitialized = .true.
  end subroutine UTPipeline_initComm


  subroutine UTPipeline_finalizeComm(doAsyncReturn)
    implicit none
    logical, optional, intent(IN) :: doAsyncReturn
    integer :: i, ierr
    logical :: doSyncReturn

    if (utpipe_numChannels > 0) then
       do i = 1, utpipe_numChannels
          if (utpipe_recvRequest(i) /= MPI_REQUEST_NULL) then
             call MPI_Cancel(utpipe_recvRequest(i), ierr)
             ASSERT_MPI_SUCCESS(ierr)
          end if
       end do
       call MPI_Waitall(utpipe_numChannels, utpipe_recvRequest, &
            utpipe_recvStatus, ierr)
       ASSERT_MPI_SUCCESS(ierr)
       utpipe_isRecvCommDone = .true.


       do i = 1, utpipe_numChannels
          if (utpipe_sendRequest(i) /= MPI_REQUEST_NULL) then
             call MPI_Cancel(utpipe_sendRequest(i), ierr)
             ASSERT_MPI_SUCCESS(ierr)
          end if
       end do
       call MPI_Waitall(utpipe_numChannels, utpipe_sendRequest, &
            utpipe_sendStatus, ierr)
       ASSERT_MPI_SUCCESS(ierr)
       utpipe_sendState(:) = CLOSE_STATE
       utpipe_isSendCommDone = .true.
    end if

    if (utpipe_size > 1) then
       if (present(doAsyncReturn)) then
          doSyncReturn = .not.doAsyncReturn
       else
          doSyncReturn = .true.
       end if
       if (doSyncReturn) then
          call MPI_Barrier(utpipe_comm, ierr)
          ASSERT_MPI_SUCCESS(ierr)
       end if
    end if
    utpipe_isCommDone = .true.
    utpipe_isCommInitialized = .false.
  end subroutine UTPipeline_finalizeComm


  subroutine UTPipeline_finalize
    implicit none
    if (utpipe_isCommInitialized) then
       call UTPipeline_finalizeComm(doAsyncReturn=.true.)
    end if

    if (utpipe_isInitialized .and. utpipe_numChannels > 0) then
       deallocate(utpipe_itemBuf)
       deallocate(utpipe_sendBuf)
       deallocate(utpipe_sendStatus)
       deallocate(utpipe_sendRequest)
       deallocate(utpipe_sendIndex)
       deallocate(utpipe_sendCount)
       deallocate(utpipe_sendState)
       deallocate(utpipe_recvBuf)
       deallocate(utpipe_recvStatus)
       deallocate(utpipe_recvRequest)
       deallocate(utpipe_recvIndex)
       deallocate(utpipe_recvCount)
       deallocate(utpipe_procList)
       utpipe_isInitialized = .false.
    end if
  end subroutine UTPipeline_finalize


  subroutine UTPipeline_progressRecvComm
    implicit none
    integer :: outcount, index, ierr, i, msgLen, procID
    logical :: isSaved

    if (utpipe_numChannels > 0 .and. .not.utpipe_isRecvCommDone) then
       !Copy old items from receive channels in which there are no pending MPI
       !receives.  This only happens when there was previously insufficient
       !space in utpipe_itemBuf.
       call utpipe_handleOldRecvMsg()

       !Test all receive channels for new messages.  Save the corresponding
       !items and then post a new receive.
       call MPI_Testsome(utpipe_NumChannels, utpipe_recvRequest, outcount, &
            utpipe_recvIndex, utpipe_recvStatus, ierr)
       ASSERT_MPI_SUCCESS(ierr)

       do i = 1, outcount
          index = utpipe_recvIndex(i)
          procID = utpipe_recvStatus(MPI_SOURCE,i)
          if (procID /= utpipe_procList(index)) then
             call Driver_abortFlash("ProcID mismatch")
          end if

          call MPI_Get_count(utpipe_recvStatus(:,i), FLASH_REAL, msgLen, ierr)
          ASSERT_MPI_SUCCESS(ierr)
          utpipe_recvCount(index) = msgLen / utpipe_itemSize
          if (utpipe_doLog) then
             write(utpipe_logUnit,'(2(a,i6))') 'Received ', &
                  utpipe_recvCount(index), ' elements from ', procID
          end if

          !A zero-byte message means the send channel is closed.
          if (msgLen > 0) then
             call utpipe_saveRecvItems(index, isSaved)
             if (isSaved) call utpipe_postRecvMsg(index)
          end if
       end do

       !Check for completion
       utpipe_isRecvCommDone = &
            all(utpipe_recvRequest == MPI_REQUEST_NULL .and. &
            utpipe_recvCount == 0)
    end if
  end subroutine UTPipeline_progressRecvComm


  subroutine utpipe_handleOldRecvMsg()
    implicit none
    integer :: i
    logical :: isSaved
    do i = 1, utpipe_numChannels
       if ( utpipe_recvRequest(i) == MPI_REQUEST_NULL .and. &
            utpipe_recvCount(i) > 0 ) then
          call utpipe_saveRecvItems(i, isSaved)
          if (isSaved) call utpipe_postRecvMsg(i)
       end if
    end do
  end subroutine utpipe_handleOldRecvMsg


  subroutine UTPipeline_sendFullestChannel
    implicit none
    integer, parameter :: notFound = -1 !This value must be negative
    integer :: fullestChannel, bufSize, i

    if (utpipe_numChannels > 0) then
       if (any(utpipe_sendCount(:) > 0)) then
          fullestChannel = notFound
          bufSize = notFound
          do i = 1, utpipe_numChannels
             !Test for data that is not currently being sent
             if ( utpipe_sendState(i) == OPEN_STATE .and. &
                  utpipe_sendRequest(i) == MPI_REQUEST_NULL .and. &
                  utpipe_sendCount(i) > 0 ) then
                if (utpipe_sendCount(i) > bufSize) then
                   fullestChannel = i
                   bufSize = utpipe_sendCount(i)
                end if
             end if
          end do
          !It is possible that there are no sends meeting the above criteria
          if (fullestChannel >= 1 .and. fullestChannel <= utpipe_numChannels) then
             call utpipe_postSendMsg(fullestChannel)
          end if
       end if
    end if
  end subroutine UTPipeline_sendFullestChannel


  subroutine utpipe_postSendMsg(index)
    implicit none
    integer, intent(IN) :: index
    integer :: procID, msgSize, ierr

    msgSize = utpipe_sendCount(index)
    if (msgSize >= 0) then
       procID = utpipe_procList(index)
       if (utpipe_doLog) then
          write(utpipe_logUnit,'(2(a,i6))') 'Post send msg of ', &
               msgSize, ' to ', procID
       end if
       call MPI_Isend(utpipe_sendBuf(1,1,index), utpipe_itemSize*msgSize, &
            FLASH_REAL, procID, utpipe_tag, utpipe_comm, &
            utpipe_sendRequest(index), ierr)
       ASSERT_MPI_SUCCESS(ierr)
    end if
    utpipe_isSendCommDone = .false.
  end subroutine utpipe_postSendMsg


  subroutine utpipe_postRecvMsg(index)
    implicit none
    integer, intent(IN) :: index
    integer :: procID, ierr

    procID = utpipe_procList(index)
    if (utpipe_doLog) then
       write(utpipe_logUnit,'(a,i6)') 'Post receive msg from ', procID
    end if
    call MPI_Irecv(utpipe_recvBuf(1,1,index), &
         utpipe_itemSize*utpipe_channelSize, FLASH_REAL, procID, &
         utpipe_tag, utpipe_comm, utpipe_recvRequest(index), ierr)
    ASSERT_MPI_SUCCESS(ierr)
    utpipe_isRecvCommDone = .false.
  end subroutine utpipe_postRecvMsg


  subroutine UTPipeline_progressComm(doFlush)
    implicit none
    logical, optional, intent(IN) :: doFlush
    call UTPipeline_progressRecvComm()
    call UTPipeline_progressSendComm()

    !The following code guarantees global progress.  It will normally
    !be executed when we are processing the last few items
    if (present(doFlush)) then
       if (doFlush .and. utpipe_itemCount == 0) then
          call UTPipeline_sendFullestChannel()
       end if
    end if
  end subroutine UTPipeline_progressComm


  !Rename UTPipeline_progressSendComm
  subroutine UTPipeline_progressSendComm()
    implicit none
    integer :: outcount, index, ierr, i

    if (utpipe_numChannels > 0 .and. .not.utpipe_isSendCommDone) then

       call utpipe_progressClosePromise()

       call MPI_Testsome(utpipe_NumChannels, utpipe_sendRequest, &
            outcount, utpipe_sendIndex, utpipe_sendStatus, ierr)
       ASSERT_MPI_SUCCESS(ierr)

       !Note that status objects are only meaningful for receive messages.
       do i = 1, outcount
          index = utpipe_sendIndex(i)

          utpipe_sendCount(index) = 0
          if (utpipe_sendState(index) == WAITING_TO_CLOSE_STATE) then
             utpipe_sendState(index) = CLOSE_STATE
          end if
          if (utpipe_doLog) then
             write(utpipe_logUnit,'(a,i6)') 'Completed send msg to ', &
                  utpipe_procList(index)
          end if
       end do

       !Check for completion
       utpipe_isSendCommDone = all(utpipe_sendState == CLOSE_STATE)
       if (utpipe_isSendCommDone .and. &
            any(utpipe_sendRequest /= MPI_REQUEST_NULL .or. &
            utpipe_sendCount /= 0)) then
          call Driver_abortFlash('Bad shutdown')
       end if
    end if
  end subroutine UTPipeline_progressSendComm


  !Promise to close the send channels.
  subroutine UTPipeline_closeSendChannels(isClosing)
    implicit none
    logical, intent(OUT) :: isClosing
    integer :: i
    if (utpipe_numChannels > 0) then
       do i = 1, utpipe_numChannels
          if (utpipe_sendState(i) == OPEN_STATE) then
             utpipe_sendState(i) = PROMISE_TO_CLOSE_STATE
          end if
       end do
       call UTPipeline_progressSendComm()
    end if
    isClosing = .true.
  end subroutine UTPipeline_closeSendChannels


  !Fulfill the close promise by sending a zero-byte notification message
  subroutine utpipe_progressClosePromise()
    implicit none
    integer :: i
    do i = 1, utpipe_numChannels
       if ( utpipe_sendState(i) == PROMISE_TO_CLOSE_STATE .and. &
            utpipe_sendRequest(i) == MPI_REQUEST_NULL ) then
          utpipe_sendCount(i) = 0 !For a zero-byte message
          call utpipe_postSendMsg(i)
          utpipe_sendState(i) = WAITING_TO_CLOSE_STATE
       end if
    end do
  end subroutine utpipe_progressClosePromise


  !Call this after UTPipeline_closeSendChannels
  subroutine UTPipeline_isCommDone(isCommDone)
    implicit none
    logical, intent(OUT) :: isCommDone

    if (.not.utpipe_isCommDone) then
       if (utpipe_isSendCommDone .and. utpipe_isRecvCommDone) then
          call UTPipeline_finalizeComm(doAsyncReturn=.true.)
       else
          call UTPipeline_progressComm()
       end if
    end if
    isCommDone = utpipe_isCommDone
  end subroutine UTPipeline_isCommDone


  subroutine UTPipeline_isDone(isDone)
    implicit none
    logical, intent(OUT) :: isDone
    logical :: isCommDone

    call UTPipeline_isCommDone(isCommDone)
    isDone = isCommDone .and. utpipe_itemCount == 0
  end subroutine UTPipeline_isDone


  subroutine UTPipeline_numItems(numItems)
    implicit none
    integer, intent(OUT) :: numItems
    numItems = utpipe_itemCount
  end subroutine UTPipeline_numItems


  !Remove this subroutine
  subroutine UTPipeline_isEmpty(isEmpty)
    implicit none
    logical, intent(OUT) :: isEmpty
    logical :: isCommBufEmpty

    if (utpipe_numChannels > 0) then
       isCommBufEmpty = &
            all(utpipe_sendCount(:) == 0 .and. utpipe_recvCount(:) == 0)
    else
       isCommBufEmpty = .true.
    end if
    isEmpty = (utpipe_itemCount == 0 .and. isCommBufEmpty)
  end subroutine UTPipeline_isEmpty


  subroutine UTPipeline_getItems(userArray, userMaxCount, userCount)
    implicit none
    real, dimension(:,:), intent(INOUT) :: userArray
    integer, intent(IN) :: userMaxCount
    integer, intent(INOUT) :: userCount
    integer :: freeSpace, itemsToCopy, firstItem

    freeSpace = userMaxCount - userCount
    itemsToCopy = min(utpipe_itemCount, freeSpace)
    if (itemsToCopy > 0) then
       !Copy from the end of utpipe_itemBuf to allow for a fast memcpy
       firstItem = utpipe_itemCount - itemsToCopy + 1

       if (utpipe_doLog) then
          write(utpipe_logUnit,'(2(a,2(i6)))') 'Copy from buf slice ', &
               firstItem, firstItem+itemsToCopy-1, ' to user slice ', &
               userCount+1, userCount+itemsToCopy
       end if

       userArray(:,userCount+1:userCount+itemsToCopy) = &
            utpipe_itemBuf(:,firstItem:firstItem+itemsToCopy-1)

       utpipe_itemCount = utpipe_itemCount - itemsToCopy
       userCount = userCount + itemsToCopy

       if (utpipe_doLog) then
          write(utpipe_logUnit,'(2(a,i6))') 'Elements in user array ', &
               userCount, ' Elements in buf ', utpipe_itemCount
       end if
    end if
  end subroutine UTPipeline_getItems


  subroutine utpipe_saveRecvItems(index, isSaved)
    implicit none
    integer, intent(IN) :: index
    logical, intent(OUT) :: isSaved
    integer :: numItems

    numItems = utpipe_recvCount(index)
    if (numItems > 0) then
       if (utpipe_itemCount + numItems <= utpipe_maxItems) then

          if (utpipe_doLog) then
             write(utpipe_logUnit,'(a,i6,2(a,2(i6)))') 'Handle receive from ', &
                  utpipe_procList(index), ' copy from msg slice ', &
                  1, numItems, ' to buf slice ', utpipe_itemCount+1, &
                  utpipe_itemCount+numItems
          end if

          utpipe_itemBuf(:,utpipe_itemCount+1:utpipe_itemCount+numItems) = &
               utpipe_recvBuf(:,1:numItems,index)
          utpipe_itemCount = utpipe_itemCount + numItems
          utpipe_recvCount(index) = 0
          isSaved = .true.
       else
          isSaved = .false.
       end if
    end if
  end subroutine utpipe_saveRecvItems


  !Caller should probably add the following: if (isHandled) item = NONEXISTENT
  subroutine UTPipeline_sendItem(item, procID, isHandled)
    implicit none
    real, dimension(:), intent(IN) :: item
    integer, intent(IN) :: procID
    logical, intent(OUT) :: isHandled
    integer :: channel, ptr, i
    integer, parameter :: notFound = -1

    !It may be necessary to change the utpipe_procList data structure
    !to make the lookup faster.
    channel = notFound
    do i = 1, utpipe_numChannels
       if (utpipe_procList(i) == procID) then
          channel = i
          exit
       end if
    end do
    if (channel == notFound) call Driver_abortFlash("Msg channel not found")

    !If there is a pending send in our desired channel we test all
    !send channels.  Request values are reset to MPI_REQUEST_NULL when
    !sends complete.
    if (utpipe_sendRequest(channel) /= MPI_REQUEST_NULL) then
       call UTPipeline_progressSendComm()
    end if

    !We can safetly add items to the send buffer if there is no pending send.
    if ( utpipe_sendState(channel) == OPEN_STATE .and. &
         utpipe_sendRequest(channel) == MPI_REQUEST_NULL ) then
       ptr = utpipe_sendCount(channel) + 1
       if (ptr > utpipe_channelSize) call Driver_abortFlash("Counting error")
       utpipe_sendBuf(:,ptr,channel) = item(:)
       utpipe_sendCount(channel) = ptr !Array is needed in utpipe_postSendMsg
          
       if (utpipe_sendCount(channel) == utpipe_channelSize) then
          call utpipe_postSendMsg(channel)
       end if
       isHandled = .true.
    else
       isHandled = .false.
    end if
  end subroutine UTPipeline_sendItem


  !We allow the caller to see the internal state of the pipeline
  !message exchange.
  subroutine UTPipeline_iterateItems(readOnlyFn)
    implicit none
    interface
       subroutine readOnlyFn(item, itemDescription)
         implicit none
         real, dimension(:), intent(IN) :: item
         character(len=*), intent(IN) :: itemDescription
       end subroutine readOnlyFn
    end interface
    integer :: i, n
    character(len=100) :: itemDescription

    do i = 1, utpipe_itemCount
       write (itemDescription,'(a,i10)') 'itemBuf: item ', i
       call readOnlyFn(utpipe_itemBuf(:,i), trim(itemDescription))
    end do

    !"The sender should not modify any part of the send buffer after a
    !nonblocking send operation is called, until the send completes."
    ![MPI-3 3.7.2].  "(the send operation itself leaves the content of
    !the send buffer unchanged)" [MPI-3 3.7.3]
    !... it should therefore always be OK to read what is there.
    do n = 1, utpipe_numChannels
       if (utpipe_sendCount(n) > 0) then
          do i = 1, utpipe_sendCount(n)
             write (itemDescription,'(2(a,i10))') 'sendBuf: channel ', n, &
                  & ', item ', i
             call readOnlyFn(utpipe_sendBuf(:,i,n), trim(itemDescription))
          end do
       end if
    end do

    !"The receiver should not access any part of the receive buffer
    !after a nonblocking receive operation is called, until the
    !receive completes." [MPI-3 3.7.2].
    !... it is not OK to read what is there.
    do n = 1, utpipe_numChannels
       if ( utpipe_recvCount(n) > 0 .and. &
            utpipe_recvRequest(n) == MPI_REQUEST_NULL ) then
          do i = 1, utpipe_recvCount(n)
             write (itemDescription,'(2(a,i10))') 'recvBuf: channel ', n, &
                  & ', item ', i
             call readOnlyFn(utpipe_recvBuf(:,i,n), trim(itemDescription))
          end do
       end if
    end do
  end subroutine UTPipeline_iterateItems


#ifdef UTPIPELINE_UNIT_TEST
  subroutine Driver_checkMPIErrorCode(errorCode)
    implicit none
    integer, intent(IN) :: errorCode
    if (errorCode /= MPI_SUCCESS) call Driver_abortFlash('Error in MPI')
  end subroutine Driver_checkMPIErrorCode

  subroutine Driver_abortFlash(msg)
    implicit none
    character (len=*), intent(IN) :: msg
    print *, "ERROR!!! ", msg
    stop
  end subroutine Driver_abortFlash
#endif

end module UTPipeline
