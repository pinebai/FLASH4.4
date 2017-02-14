!!****ih* source/flashUtilities/UTCounter/MasterSlave/UTCounter_sharedCounter
!!
!!  NAME
!!    UTCounter_sharedCounter
!!
!!  SYNOPSIS
!!    UTCounter_sharedCounter is a fortran module
!!
!!  DESCRIPTION
!!    Utility module which maintains a single shared counter over all MPI ranks
!!    in a passed MPI communicator.  All MPI ranks are able to increment
!!    the counter and will be notified when the counter target is met.
!!
!!  ARGUMENTS
!!    none
!!
!!  EXAMPLE
!!    none
!!
!!  NOTES
!!    The shared counter is limited by the maximum value of a default integer
!!    which is generally 2,147,483,647.  If the count is to be greater than
!!    the maximum value of a default integer then the portable way to change
!!    the code is to introduce a new integer type named "utcnt_int":
!!
!!    module UTCounter_sharedCounterType
!!      implicit none
!!      !Minimum number of digits in the integer type.  The table is generally:
!!      ! Digits Bytes
!!      ! 2      1
!!      ! 4      2
!!      ! 9      4     <-- This is normally the default
!!      ! 18     8
!!      ! 38     16
!!      integer, parameter :: utcnt_int_digits = 18
!!      !The counter integer type which can hold the required number of digits
!!      integer, parameter :: utcnt_int = selected_int_kind(utcnt_int_digits)
!!    end module UTCounter_sharedCounterType
!!
!!    We can access the integer and appropriate MPI type from anywhere:
!!    use UTCounter_sharedCounterType
!!    integer(kind=utcnt_int) :: counterTarget
!!    integer :: utcnt_int_mpi_type
!!    call MPI_Type_create_f90_integer(utcnt_int_digits, utcnt_int_mpi_type, ierr)
!!
!!***

module UTCounter_sharedCounter

#ifdef UTCOUNTER_UNIT_TEST
  use mpi
#else
  use Driver_interface, ONLY : Driver_abortFlash, Driver_checkMPIErrorCode
#endif

  implicit none
  private
  public :: UTCounter_init, UTCounter_startCounter, UTCounter_incrementCounter, &
       UTCounter_progressCounter, UTCounter_stopCounter, UTCounter_finalize

#ifdef UTCOUNTER_UNIT_TEST
  integer, parameter :: FLASH_INTEGER = MPI_INTEGER
#else
  include 'Flash_mpi.h'
#endif

  integer, save :: utcnt_counterTarget
  integer, save :: utcnt_sharedCounter
  integer, save :: utcnt_localCounter

  integer, save :: utcnt_comm
  integer, save :: utcnt_rank
  integer, save :: utcnt_size
  integer, save :: utcnt_masterRank
  
  integer, parameter :: utcnt_contributorTag = 1236
  integer, allocatable, save, dimension(:,:) :: utcnt_contributorStatus
  integer, allocatable, save, dimension(:) :: utcnt_contributorRequest
  integer, allocatable, save, dimension(:) :: utcnt_contributorIndex
  integer, allocatable, save, dimension(:) :: utcnt_contributorMsg
  integer, save :: utcnt_contributorSize

  integer, parameter :: utcnt_sharedCountTag = 1237
  integer, allocatable, save, dimension(:,:) :: utcnt_sharedCountStatus
  integer, allocatable, save, dimension(:) :: utcnt_sharedCountRequest
  integer, save :: utcnt_sharedCountSize
  integer, save :: utcnt_sharedCountMsg

  integer, save :: utcnt_logUnit
  logical, save :: utcnt_doLog

  logical, save :: utcnt_isInitialized = .false.
  logical, save :: utcnt_isCounterTargetMet = .false.

  !utcnt_isDone is .true. when all outstanding communication is complete.
  logical, save :: utcnt_isDone = .true.

#ifdef FLASH_MPI3
  logical, parameter :: utcnt_useMPI3 = .true.
#else
  logical, parameter :: utcnt_useMPI3 = .false.
#endif

#ifdef FLASH_LIBNBC
# define NBC_OK 0
# define NBC_CONTINUE 3
# define NBC_REQUEST_NULL -1
# define FLASH_REQUEST_NULL NBC_REQUEST_NULL
#else
# define FLASH_REQUEST_NULL MPI_REQUEST_NULL
#endif

  integer, parameter :: utcnt_invalid_count = -1

contains


  !!****if* source/flashUtilities/UTCounter/UTCounter_init
  !!
  !!  NAME     
  !!   UTCounter_init
  !!
  !!  SYNOPSIS
  !!   UTCounter_init(integer(IN)           :: comm,
  !!                  optional, integer(IN) :: masterRank,
  !!                  optional, integer(IN) :: logUnit)
  !!
  !!  DESCRIPTION 
  !!    Initializes the counter module.  All MPI ranks must pass the same
  !!    communicator and master rank in the subroutine call.
  !!
  !!  ARGUMENTS
  !!    comm: MPI communicator in which we will maintain the counter
  !!    masterRank: The master MPI rank
  !!    logUnit: An optional Fortran unit number for writing log messages
  !!
  !!***

  subroutine UTCounter_init(comm, masterRank, logUnit)
    implicit none
    integer, intent(IN) :: comm
    integer, optional, intent(IN) :: masterRank, logUnit
    integer :: ierr

    utcnt_comm = comm
    call MPI_Comm_rank(utcnt_comm, utcnt_rank, ierr)
    call Driver_checkMPIErrorCode(ierr)
    call MPI_Comm_size(utcnt_comm, utcnt_size, ierr)
    call Driver_checkMPIErrorCode(ierr)

    utcnt_masterRank = 0
    if (present(masterRank)) then
       utcnt_masterRank = masterRank
    end if
    if (utcnt_masterRank < 0 .or. utcnt_masterRank >= utcnt_size) then
       call Driver_abortFlash('Invalid Master')
    end if

    utcnt_logUnit = -1
    utcnt_doLog = .false.
    if (present(logUnit)) then
       utcnt_doLog = (logUnit >= 0) !Must be a valid Fortran unit
       if (utcnt_doLog) utcnt_logUnit = logUnit
    end if

    if (utcnt_size > 1) then
       if (utcnt_rank == utcnt_masterRank) then
          utcnt_contributorSize = utcnt_size - 1
          if (utcnt_useMPI3) then
             utcnt_sharedCountSize = 1
          else
             utcnt_sharedCountSize = utcnt_size - 1
          endif
       else
          utcnt_contributorSize = 1
          utcnt_sharedCountSize = 1
       end if

       allocate(utcnt_contributorStatus(MPI_STATUS_SIZE,utcnt_contributorSize))
       allocate(utcnt_contributorRequest(utcnt_contributorSize))
       allocate(utcnt_contributorIndex(utcnt_contributorSize))
       allocate(utcnt_contributorMsg(utcnt_contributorSize))
       allocate(utcnt_sharedCountStatus(MPI_STATUS_SIZE,utcnt_sharedCountSize))
       allocate(utcnt_sharedCountRequest(utcnt_sharedCountSize))
    end if
    utcnt_isInitialized = .true.
  end subroutine UTCounter_init


  !!****if* source/flashUtilities/UTCounter/UTCounter_startCounter
  !!
  !!  NAME     
  !!   UTCounter_startCounter
  !!
  !!  SYNOPSIS
  !!   UTCounter_startCounter(integer(IN) :: counterTarget)
  !!
  !!  DESCRIPTION 
  !!    Initializes the communication which will persist until either
  !!      1) the counter target is met, or
  !!      2) the subroutine UTCounter_stopCounter is called
  !!
  !!  ARGUMENTS
  !!    counterTarget: The actual counter target
  !!
  !!***

  subroutine UTCounter_startCounter(counterTarget)
    implicit none
    integer, intent(IN) :: counterTarget
    integer :: index, procID, ierr

    if (.not.utcnt_isInitialized) then
       call Driver_abortFlash("Must intialize the counter first")
    end if

    utcnt_counterTarget = counterTarget
    utcnt_localCounter = 0
    utcnt_sharedCounter = 0

    utcnt_isCounterTargetMet = utcnt_sharedCounter >= utcnt_counterTarget
    if (utcnt_isCounterTargetMet) then
       utcnt_sharedCountMsg = utcnt_sharedCounter
    else
       utcnt_sharedCountMsg = utcnt_invalid_count
    end if

    utcnt_isDone = utcnt_isCounterTargetMet
    if (utcnt_size > 1 .and. .not.utcnt_isDone) then
       if (utcnt_rank == utcnt_masterRank) then
          index = 1
          do procID = 0, utcnt_size-1
             if (procID /= utcnt_masterRank) then
                call MPI_Irecv(utcnt_contributorMsg(index), 1, FLASH_INTEGER, &
                        procID, utcnt_contributorTag, utcnt_comm, &
                        utcnt_contributorRequest(index), ierr)
                call Driver_checkMPIErrorCode(ierr)
                index = index + 1
             end if
          end do
       else
          !Must be initialized to MPI_REQUEST_NULL because slaves check whether
          !there is an active contribution before they send a new contribution.
          utcnt_contributorRequest(1) = MPI_REQUEST_NULL
          if (utcnt_useMPI3) then
             call utcnt_ibcastDone()
          else
             call MPI_Irecv(utcnt_sharedCountMsg, 1, FLASH_INTEGER, &
                  utcnt_masterRank, utcnt_sharedCountTag, utcnt_comm, &
                  utcnt_sharedCountRequest, ierr)
             call Driver_checkMPIErrorCode(ierr)
          endif
       end if
    end if
  end subroutine UTCounter_startCounter


  subroutine utcnt_ibcastDone()
    implicit none
#ifdef FLASH_MPI3
    integer :: ierr
# ifdef FLASH_LIBNBC
    !LibNBC NOTE: Fortran types are *not* supported.  We assume
    !sizeof(MPI_INT) == sizeof(MPI_INTEGER) == sizeof(FLASH_INTEGER).
    call NBC_Ibcast(utcnt_sharedCountMsg, 1, MPI_INT, utcnt_masterRank, &
         utcnt_comm, utcnt_sharedCountRequest, ierr)
    if (ierr /= NBC_OK) call Driver_abortFlash("NBC_Ibcast error")
# else
    call MPI_Ibcast(utcnt_sharedCountMsg, 1, FLASH_INTEGER, utcnt_masterRank, &
         utcnt_comm, utcnt_sharedCountRequest, ierr)
    call Driver_checkMPIErrorCode(ierr)
# endif
#endif
  end subroutine utcnt_ibcastDone


  subroutine utcnt_waitDone
    implicit none
    integer :: ierr

#ifdef FLASH_LIBNBC
    !LibNBC NOTE: The Fortran versions of NBC_Test / NBC_Wait do *not*
    !initialize the error code for NBC_REQUEST_NULL requests.
    if (utcnt_sharedCountRequest(1) /= FLASH_REQUEST_NULL) then
       call NBC_Wait(utcnt_sharedCountRequest(1), ierr)
       if (ierr /= NBC_OK) call Driver_abortFlash("NBC_Wait error")
    end if
#else
    call MPI_Wait(utcnt_sharedCountRequest(1), &
         utcnt_sharedCountStatus(:,1), ierr)
    call Driver_checkMPIErrorCode(ierr)
#endif
  end subroutine utcnt_waitDone


  subroutine utcnt_testDone(isDone)
    implicit none
    logical, intent(OUT) :: isDone
    integer :: ierr

    isDone = .false.
    if (utcnt_sharedCountRequest(1) == FLASH_REQUEST_NULL) then
       call Driver_abortFlash("Should not be null - test")
    end if
#ifdef FLASH_LIBNBC
    !LibNBC NOTE: The Fortran versions of NBC_Test / NBC_Wait do *not*
    !initialize the error code for NBC_REQUEST_NULL requests.
    call NBC_Test(utcnt_sharedCountRequest(1), ierr)
    if (ierr /= NBC_OK .and. ierr /= NBC_CONTINUE) &
         call Driver_abortFlash("NBC_Test error")
    isDone = (utcnt_sharedCountRequest(1) == FLASH_REQUEST_NULL)
#else
    call MPI_Test(utcnt_sharedCountRequest(1), isDone, &
         utcnt_sharedCountStatus(:,1), ierr)
    call Driver_checkMPIErrorCode(ierr)
#endif
  end subroutine utcnt_testDone



  !isTargetMet is only ever set to .true. once all outstanding
  !communications are either complete or cancelled, that is
  !utcnt_isDone is .true..  It is impossible for isTargetMet to be
  !.true. and utcnt_isDone to be .false.

  subroutine UTCounter_progressCounter(isTargetMet, finalCount)
    implicit none
    logical, intent(OUT) :: isTargetMet
    integer, optional, intent(OUT) :: finalCount

    if (utcnt_rank == utcnt_masterRank) then
       call utcnt_progressMasterCounter(isTargetMet)
    else
       call utcnt_progressSlaveCounter(isTargetMet)
    end if

    if (isTargetMet .and. .not.utcnt_isDone) &
         call Driver_abortFlash("Progress error")

    if (present(finalCount)) then
       if (isTargetMet) then
          !It only makes sense to return the final count *after* the
          !target is met.  Do not access utcnt_sharedCountMsg *before*
          !the target is met because there will be an outstanding
          !MPI_Ibcast or MPI_Isend / MPI_Irecv on utcnt_sharedCountMsg.
          finalCount = utcnt_sharedCountMsg
       else
          finalCount = utcnt_invalid_count
       end if
    end if
  end subroutine UTCounter_progressCounter


  subroutine utcnt_progressMasterCounter(isTargetMet)
    implicit none
    logical, intent(OUT) :: isTargetMet
    integer :: outcount, procID, index, ierr, i

    if (.not.utcnt_isDone) then
       !Add contribution from myself (the master)
       if (utcnt_localCounter > 0) then
          call utcnt_logContribution(utcnt_rank, utcnt_masterRank, &
               utcnt_localCounter)
          utcnt_sharedCounter = utcnt_sharedCounter + utcnt_localCounter
          utcnt_localCounter = 0
       end if

       !Add contribution from slaves
       if (utcnt_size > 1 .and. utcnt_sharedCounter < utcnt_counterTarget) then
          call MPI_Testsome(utcnt_size-1, utcnt_contributorRequest, &
               outcount, utcnt_contributorIndex, utcnt_contributorStatus, &
               ierr)
          call Driver_checkMPIErrorCode(ierr)

          if (outcount >= 1) then
             do i = 1, outcount
                procID = utcnt_contributorStatus(MPI_SOURCE,i)
                index = utcnt_contributorIndex(i)
                call utcnt_logContribution(procID, utcnt_rank, &
                     utcnt_contributorMsg(index))

                utcnt_sharedCounter = utcnt_sharedCounter + &
                     utcnt_contributorMsg(index)
                utcnt_contributorMsg(index) = 0

                !Post a new speculative receive
                call MPI_Irecv(utcnt_contributorMsg(index), 1, FLASH_INTEGER, &
                     procID, utcnt_contributorTag, utcnt_comm, &
                     utcnt_contributorRequest(index), ierr)
                call Driver_checkMPIErrorCode(ierr)
             end do
          end if
       end if

       !Check if we have met the counter target.  If so, inform the slaves
       utcnt_isCounterTargetMet = utcnt_sharedCounter >= utcnt_counterTarget
       if (utcnt_isCounterTargetMet) then
          utcnt_sharedCountMsg = utcnt_sharedCounter
          if (utcnt_size > 1) then
             if (utcnt_useMPI3) then       
                call utcnt_ibcastDone()
                call utcnt_waitDone()
             else
                index = 1
                do procID = 0, utcnt_size-1
                   if (procID /= utcnt_masterRank) then
                      call MPI_Isend(utcnt_sharedCountMsg, 1, FLASH_INTEGER, &
                           procID, utcnt_sharedCountTag, utcnt_comm, &
                           utcnt_sharedCountRequest(index), ierr)
                      call Driver_checkMPIErrorCode(ierr)
                      index = index + 1
                   end if
                end do
                call MPI_Waitall(utcnt_size-1, utcnt_sharedCountRequest, &
                     utcnt_sharedCountStatus, ierr)
                call Driver_checkMPIErrorCode(ierr)
             end if
          end if
          !utcnt_isDone is set to .true. in UTCounter_stopCounter
          call UTCounter_stopCounter(doAsyncReturn=.true.)
       end if
    end if !END: if (.not.utcnt_isDone)
    isTargetMet = utcnt_isCounterTargetMet
  end subroutine utcnt_progressMasterCounter


  subroutine utcnt_progressSlaveCounter(isTargetMet)
    implicit none
    logical, intent(OUT) :: isTargetMet
    integer :: ierr
    logical :: flag

    if (.not.utcnt_isDone) then
       if (utcnt_localCounter > 0) then
          !The MPI_REQUEST_NULL initialization is needed so that it is safe
          !to call MPI_Test before we have actually sent a message:
          !"One is allowed to call MPI_TEST with a null or inactive request
          !argument. In such a case the operation returns with flag = true
          !and empty status." - MPI-3.0 standard, section 3.7.3.
          call MPI_Test(utcnt_contributorRequest(1), flag, &
               utcnt_contributorStatus(:,1), ierr)
          call Driver_checkMPIErrorCode(ierr)

          if (flag) then
             call utcnt_logContribution(utcnt_rank, utcnt_masterRank, &
                  utcnt_localCounter)
             utcnt_contributorMsg(1) = utcnt_localCounter
             utcnt_localCounter = 0

             !I use MPI_Issend instead of MPI_Isend to reduce the number of
             !messages that each slave sends to the master: "[An MPI_Wait on
             !an MPI_Isend]... does not indicate that the message has been
             !received, rather, it may have been buffered by the communication
             !subsystem. However, if a synchronous mode send was used, the
             !completion of the send operation indicates that a matching
             !receive was initiated, and that the message will eventually be
             !received by this matching receive" - MPI-3.0 standard,
             !section 3.7.3.
             call MPI_Issend(utcnt_contributorMsg, 1, FLASH_INTEGER, &
                  utcnt_masterRank, utcnt_contributorTag, utcnt_comm, &
                  utcnt_contributorRequest, ierr)
             call Driver_checkMPIErrorCode(ierr)
          end if
       end if

       if (utcnt_useMPI3) then
          call utcnt_testDone(flag)
       else
          call MPI_Test(utcnt_sharedCountRequest(1), flag, &
               utcnt_sharedCountStatus(:,1), ierr)
          call Driver_checkMPIErrorCode(ierr)
       end if

       if (flag) then
          !utcnt_isDone is set to .true. in UTCounter_stopCounter
          call UTCounter_stopCounter(doAsyncReturn=.true.)
          utcnt_isCounterTargetMet = &
               (utcnt_sharedCountMsg /= utcnt_invalid_count)
       end if
    end if !END: if (.not.utcnt_isDone)
    isTargetMet = utcnt_isCounterTargetMet
  end subroutine utcnt_progressSlaveCounter


  subroutine utcnt_logContribution(fromRank, toRank, contribution)
    implicit none
    integer, intent(IN) :: fromRank, toRank, contribution
    character(len=200) :: msg1, msg2, fullMsg

    if (utcnt_doLog) then
       write(msg1,'(a,i8)') 'contribution of ', contribution

       if (utcnt_rank == fromRank .and. utcnt_rank == toRank) then
          write(msg2,'(a,i6)') ' copied from  ', utcnt_rank
       else
          if (utcnt_rank == fromRank) then
             write(msg2,'(a,i6)') ' sent to      ', toRank
          else if (utcnt_rank == toRank) then
             write(msg2,'(a,i6)') ' received from', fromRank
          else
             call Driver_abortFlash('Invalid arguments')
          end if
       end if

       if (utcnt_logUnit == 6) then
          write(fullMsg,'(a,i6,a,a,a)') '[Rank ', utcnt_rank, ']: ', &
               trim(msg1), trim(msg2)
       else
          fullMsg = trim(msg1)//trim(msg2)
       end if
       write(utcnt_logUnit,'(a)') trim(fullMsg)
    end if
  end subroutine utcnt_logContribution


  !Calling UTCounter_stopCounter in any manner guarantees that all
  !messages are correctly matched or cancelled.

  !Calling UTCounter_stopCounter with doAsyncReturn=.false. or
  !without any arguments gives a further guarantee that it is safe the
  !use the counter module again.  This is because the MPI_Barrier
  !synchronization prevents message interception between different
  !invocations of UTCounter.

  !We do not enforce a synchronization because the experienced user
  !may wish to delay the synchronization until a later point in the
  !code to improve global progress.
 
  subroutine UTCounter_stopCounter(doAsyncReturn)
    implicit none
    logical, optional, intent(IN) :: doAsyncReturn
    integer :: i, ierr
    logical :: doSyncReturn  

    if (utcnt_size > 1) then
       doSyncReturn = .true.
       if (present(doAsyncReturn)) then
          doSyncReturn = .not.doAsyncReturn
       end if

       if (.not.utcnt_isDone) then
          do i = 1, utcnt_contributorSize
             if (utcnt_contributorRequest(i) /= MPI_REQUEST_NULL) then
                !I found that there is an issue when MPI_REQUEST_NULL gets
                !passed to MPI_Cancel.
                call MPI_Cancel(utcnt_contributorRequest(i), ierr)
                call Driver_checkMPIErrorCode(ierr)
             end if
          end do
          call MPI_Waitall(utcnt_contributorSize, utcnt_contributorRequest, &
               utcnt_contributorStatus, ierr)
          call Driver_checkMPIErrorCode(ierr)

          !This code block handles an early shut-down, i.e. before
          !a "done" message gets sent.
          if (.not.utcnt_isCounterTargetMet) then
             if (utcnt_useMPI3) then
                if (utcnt_rank == utcnt_masterRank) then
                   !"It is erroneous to call MPI_REQUEST_FREE or MPI_CANCEL
                   !for a request associated with a nonblocking collective
                   !operation" - MPI-3.0 standard, section 5.12.
                   utcnt_sharedCountMsg = utcnt_invalid_count
                   call utcnt_ibcastDone()
                end if
                !Request may already have been set to NULL in utcnt_testDone
                call utcnt_waitDone()
             else
                if (utcnt_rank /= utcnt_masterRank) then
                   if (utcnt_sharedCountRequest(1) /= MPI_REQUEST_NULL) then
                      call MPI_Cancel(utcnt_sharedCountRequest(1), ierr)
                      call Driver_checkMPIErrorCode(ierr)
                      call MPI_Wait(utcnt_sharedCountRequest(1), &
                           utcnt_sharedCountStatus(:,1), ierr)
                      call Driver_checkMPIErrorCode(ierr)
                   end if
                end if
             end if
          end if
       end if !END: if (.not.utcnt_isDone)

       if (doSyncReturn) then
          !We have correctly matched (or cancelled) the nonblocking
          !sends and receives, however, there can be message
          !interception between different invocations of UTCounter.
          !The barrier guarantees that everyone has completed
          !communication before returning and so the next invocation
          !of UTCounter is safe.
          call MPI_Barrier(utcnt_comm, ierr)
          call Driver_checkMPIErrorCode(ierr)
       end if
    end if
    utcnt_isDone = .true.
  end subroutine UTCounter_stopCounter


  subroutine UTCounter_finalize()
    implicit none
    if (.not.utcnt_isDone) then
       call UTCounter_stopCounter(doAsyncReturn=.true.)
    end if
    if (utcnt_isInitialized) then
       if (utcnt_size > 1) then
          deallocate(utcnt_contributorStatus)
          deallocate(utcnt_contributorRequest)
          deallocate(utcnt_contributorIndex)
          deallocate(utcnt_contributorMsg)
          deallocate(utcnt_sharedCountStatus)
          deallocate(utcnt_sharedCountRequest)
       end if
       utcnt_isInitialized = .false.
       utcnt_isCounterTargetMet = .false.
    end if
  end subroutine UTCounter_finalize


  subroutine UTCounter_incrementCounter(incr)
    implicit none
    integer, intent(IN) :: incr
    utcnt_localCounter = utcnt_localCounter + incr
  end subroutine UTCounter_incrementCounter


#ifdef UTCOUNTER_UNIT_TEST
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

end module UTCounter_sharedCounter

! TODO: Need to update the following comments according to the new API
!
! It is sometimes necessary to synchronize between two uses of
! UTCounter_sharedCounter.  Example 1 was actually seen.  Example 2
! could technically happen.  Both depend on the sum of local counts
! being greater than the shared target count.
!
!
! Example 1:
!
! Rank 0 (master)                                 Rank 1 (slave)
!
! UTCounter_init                                  UTCounter_init
! UTCounter_incrementCounter                      UTCounter_incrementCounter
! UTCounter_isCounterTargetMet                    UTCounter_isCounterTargetMet
!   Target is not met, but we decide                Test for completion:
!   to terminate anyway                              MPI_Test
! UTCounter_finalize                                       ^
!                                                          |
! UTCounter_init                                           |
! UTCounter_incrementCounter                               |
! UTCounter_isCounterTargetMet                             |
!   Target is met:                                         |
!    MPI_Isend --------------------------------------------|
!    MPI_Waitall
!
!
!
! Example 2:
!
! Rank 0 (master)                                 Rank 1 (slave)
!
! UTCounter_init                                  UTCounter_init
! UTCounter_incrementCounter                      UTCounter_incrementCounter
! UTCounter_isCounterTargetMet                    UTCounter_isCounterTargetMet
!   Target is met:                                  Send contribution:
!    MPI_Isend                                       MPI_Issend
!    MPI_Waitall                                           |
!    (Assume communication buffered)                       |
! UTCounter_finalize                                       |
!   Cancel contributor communication channels:             |
!    MPI_Cancel [1]                                        |
!    MPI_Waitall                                           |
!                                                          |
! UTCounter_init                                           |
! UTCounter_incrementCounter                               |
! UTCounter_isCounterTargetMet                             |
!   MPI_Testsome <-----------------------------------------|
!
!"If a receive is marked for cancellation, then it must be the case
!that either the receive completes normally, or that the receive is
!successfully cancelled, in which case no part of the receive buffer is
!altered.  Then, any matching send has to be satisfied by another
!receive." - MPI-3.0 standard, section 3.8.4. [1]
