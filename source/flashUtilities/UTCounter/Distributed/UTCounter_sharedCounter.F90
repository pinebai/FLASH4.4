!!****ih* source/flashUtilities/UTCounter/Distributed/UTCounter_sharedCounter
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
  !Include Flash.h to obtain FLASH_LIB_NBC macro
# include "Flash.h"
  use Driver_interface, ONLY : Driver_abortFlash, Driver_checkMPIErrorCode
#endif

  implicit none
  private
  public :: UTCounter_init, UTCounter_startCounter, UTCounter_incrementCounter, &
       UTCounter_progressCounter, UTCounter_stopCounter, UTCounter_finalize

#ifdef UTCOUNTER_UNIT_TEST
  integer, parameter :: FLASH_INTEGER = MPI_INTEGER
  integer, parameter :: FLASH_LOGICAL = MPI_LOGICAL
#else
  include 'Flash_mpi.h'
#endif

  integer, save :: utcnt_localCount
  integer, save :: utcnt_targetCount

  integer, save :: utcnt_localCountMsg
  integer, save :: utcnt_sharedCountMsg
  integer, save :: utcnt_request

  integer, save :: utcnt_comm
  integer, save :: utcnt_rank
  integer, save :: utcnt_masterRank
  integer, save :: utcnt_size

  logical, save :: utcnt_isInitialized = .false.
  logical, save :: utcnt_isDone = .true.
  logical, save :: utcnt_isCounterTargetMet = .false.

  integer, save :: utcnt_logUnit
  logical, save :: utcnt_doLog

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

#ifndef FLASH_MPI3
    call Driver_abortFlash("This counter implementation requires MPI-3")
#endif

    if (utcnt_isInitialized) call UTCounter_finalize()

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

    utcnt_doLog = .false.
    utcnt_logUnit = -1
    if (present(logUnit)) then
       utcnt_doLog = (logUnit >= 0) !Must be a valid Fortran unit
       if (utcnt_doLog) utcnt_logUnit = logUnit
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
  !!      2) the subroutine UTCounter_finalizeComm is called
  !!
  !!  ARGUMENTS
  !!    counterTarget: The actual counter target
  !!
  !!***

  subroutine UTCounter_startCounter(targetCount)
    implicit none
    integer, intent(IN) :: targetCount

    if (.not.utcnt_isInitialized) then
       call Driver_abortFlash("Must intialize the counter first")
    end if
    if (.not.utcnt_isDone) then
       call Driver_abortFlash("A counter invocation is already running")
    end if

    utcnt_targetCount = targetCount
    utcnt_localCount = 0
    utcnt_request = FLASH_REQUEST_NULL

    utcnt_isCounterTargetMet = utcnt_localCount >= utcnt_targetCount
    utcnt_isDone = utcnt_isCounterTargetMet
  end subroutine UTCounter_startCounter


  subroutine UTCounter_incrementCounter(incr)
    implicit none
    integer, intent(IN) :: incr
    utcnt_localCount = utcnt_localCount + incr
  end subroutine UTCounter_incrementCounter


  subroutine UTCounter_progressCounter(isTargetMet, finalCount)
    implicit none
    logical, intent(OUT) :: isTargetMet
    integer, optional, intent(OUT) :: finalCount
    integer :: ierr
    logical :: testIfTargetMet

    if (.not.utcnt_isDone) then
       if (utcnt_size > 1) then
          testIfTargetMet = .false.
          if (utcnt_request == FLASH_REQUEST_NULL) then
             !No outstanding MPI_Iallreduce so we post an MPI_Iallreduce.
             utcnt_localCountMsg = utcnt_localCount
#ifdef FLASH_MPI3
# ifdef FLASH_LIBNBC
             !LibNBC NOTE: Fortran types are *not* supported.  We assume
             !sizeof(MPI_INT) == sizeof(MPI_INTEGER) == sizeof(FLASH_INTEGER).
             call NBC_Iallreduce(utcnt_localCountMsg, utcnt_sharedCountMsg, &
                  1, MPI_INT, MPI_SUM, utcnt_comm, utcnt_request, ierr)
             if (ierr /= NBC_OK) call Driver_abortFlash("NBC_Iallreduce error")
# else
             call MPI_Iallreduce(utcnt_localCountMsg, utcnt_sharedCountMsg, &
                  1, FLASH_INTEGER, MPI_SUM, utcnt_comm, utcnt_request, ierr)
             call Driver_checkMPIErrorCode(ierr)
# endif
#endif
          else
             !Test the outstanding MPI_Iallreduce.  It must be complete for us
             !to safely test whether the counter target is met.
#ifdef FLASH_LIBNBC
             !LibNBC NOTE: The Fortran versions of NBC_Test / NBC_Wait do *not*
             !initialize the error code for NBC_REQUEST_NULL requests.
             call NBC_Test(utcnt_request, ierr)
             if (ierr /= NBC_OK .and. ierr /= NBC_CONTINUE) &
                  call Driver_abortFlash("NBC_Test error")             
             testIfTargetMet = (utcnt_request == FLASH_REQUEST_NULL)
#else
             call MPI_Test(utcnt_request, testIfTargetMet, &
                  MPI_STATUS_IGNORE, ierr)
             call Driver_checkMPIErrorCode(ierr)
#endif
          end if
       else
          utcnt_localCountMsg = utcnt_localCount
          utcnt_sharedCountMsg = utcnt_localCount
          testIfTargetMet = .true.
       end if

       if (testIfTargetMet) then
          utcnt_isCounterTargetMet = utcnt_sharedCountMsg >= utcnt_targetCount
          call utcnt_logCount()
          if (utcnt_isCounterTargetMet) then
             call UTCounter_stopCounter(doAsyncReturn=.true.)
          end if
       end if
    end if !End if (.not.utcnt_isDone)

    isTargetMet = utcnt_isCounterTargetMet
    if (present(finalCount)) then
       if (isTargetMet) then
          !It only makes sense to return the final count *after* the
          !target is met.  Do not access utcnt_sharedCountMsg *before*
          !the target is met because there will be an outstanding
          !MPI_Iallreduce on utcnt_sharedCountMsg.
          finalCount = utcnt_sharedCountMsg
       else
          finalCount = utcnt_invalid_count
       end if
    end if
  end subroutine UTCounter_progressCounter


  subroutine utcnt_logCount()
    implicit none
    character(len=200) :: msg
    if (utcnt_doLog .and. utcnt_rank == utcnt_masterRank) then
       write(msg,'(a,i12)') 'current shared count ', utcnt_sharedCountMsg
       if (utcnt_logUnit == 6) then
          write(utcnt_logUnit,'(a,i8,a,a)') '[Rank ', utcnt_rank, ']: ', &
               trim(msg)
       else
          write(utcnt_logUnit,'(a)') trim(msg)
       end if
    end if
  end subroutine utcnt_logCount


  !This implementation of the counter does not support stopping the
  !counter before the target is met.  This is because we cannot
  !guaranteee that all MPI ranks post the same number of all reduces
  !(unless we add synchronization).

  !To illustrate the problem:
  !Scenario 1
  !                                     Rank 0                   Rank 1
  ! UTCounter_progressCounter (call 1): Post reduce              Post reduce
  ! UTCounter_progressCounter (call 2): Test reduce (complete)
  !                                     Post reduce
  ! UTCounter_stopCounter:              Wait reduce              Wait reduce
  !
  !
  !Scenario 2
  !                                     Rank 0                   Rank 1
  ! UTCounter_progressCounter (call 1): Post reduce              Post reduce
  ! UTCounter_progressCounter (call 2): Test reduce (incomplete)
  ! UTCounter_stopCounter:              Wait reduce              Wait reduce
  !
  !In scenario 1, rank 0 has posted one more reduce than rank 1 and so
  !will be waiting indefinitely.

  subroutine UTCounter_stopCounter(doAsyncReturn)
    implicit none
    logical, optional, intent(IN) :: doAsyncReturn
    integer :: ierr
    logical :: doSyncReturn

    if (utcnt_size > 1) then
       doSyncReturn = .true.
       if (present(doAsyncReturn)) then
          doSyncReturn = .not.doAsyncReturn
       end if

       if (.not.utcnt_isCounterTargetMet) then
          call Driver_abortFlash("Counter does not support early termination")
       end if

       if (doSyncReturn) then
          !This synchronization option is only here for compatibility
          !with the master-slave implementation of the counter.  The 
          !distributed counter never needs to synchronize because it does
          !not support stopping the counter before the target is met.
          call MPI_Barrier(utcnt_comm, ierr)
          call Driver_checkMPIErrorCode(ierr)
       end if
    end if
    utcnt_isDone = .true.
  end subroutine UTCounter_stopCounter


  subroutine UTCounter_finalize()
    implicit none
    if (.not.utcnt_isDone) call UTCounter_stopCounter(doAsyncReturn=.true.)
    if (utcnt_isInitialized) then
       utcnt_isInitialized = .false.
       utcnt_isCounterTargetMet = .false.
    end if
  end subroutine UTCounter_finalize


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
