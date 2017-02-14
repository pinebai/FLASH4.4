Program test_shared_counter
  call test_subroutine
End Program test_shared_counter


subroutine test_subroutine
  use mpi
  use UTCounter_sharedCounter
  implicit none
  integer, parameter :: mpiMasterRank = 0
  integer :: mpiComm, mpiRank, mpiSize
  integer :: ierr

  call MPI_Init(ierr)
  mpiComm = MPI_COMM_WORLD
  call MPI_Comm_size(mpiComm, mpiSize, ierr)
  call MPI_Comm_rank(mpiComm, mpiRank, ierr)
  call MPI_Errhandler_set(mpiComm, MPI_ERRORS_RETURN, ierr)
  
  if (mpiRank == 0) print *, "Maximum global count", HUGE(ierr)
  if (mpiRank == 0) print *, "starting test1"
  call test1()
#ifdef UNMET_TARGET_OK
  if (mpiRank == 0) print *, "starting test2"
  call test2()
#endif
  if (mpiRank == 0) print *, "starting test3"
  call test3()
#ifdef UNMET_TARGET_OK
  if (mpirank == 0) print *, "starting test4"
  call test4()
#endif
  if (mpirank == 0) print *, "starting test5"
  call test5()
  if (mpirank == 0) print *, "starting test6"
  call test6()
  if (mpirank == 0) print *, "starting test7"
  call test7()
  if (mpiRank == 0) print *, "finished all tests!!"

  call MPI_Errhandler_set(mpiComm, MPI_ERRORS_ARE_FATAL, ierr)
  call MPI_Finalize(ierr)
  
contains

  !Test a simple start-up and shut-down
  subroutine test1()
    implicit none
    call UTCounter_init(mpiComm, mpimasterRank)
    call UTCounter_finalize()
  end subroutine test1


  !Test a simple start-up and shut-down with communication
  subroutine test2()
    implicit none
    integer :: counterTarget
    counterTarget = 100

    call UTCounter_init(mpiComm, mpimasterRank)
    call UTCounter_startCounter(counterTarget)
    call UTCounter_stopCounter()
    call UTCounter_finalize()
  end subroutine test2


  !Test the counter in a MPI_COMM_SELF communicator
  subroutine test3()
    implicit none
    integer, parameter :: logUnit = 23
    integer :: counterTarget
    logical :: isCommDone, isTargetMet
    character(len=100) :: fileName

    counterTarget = 20
    isCommDone = .false.
    write(fileName,'(a,i6.6,a)') 'test3_', mpirank, '.out'
    open (logUnit, file=fileName)

    call UTCounter_init(MPI_COMM_SELF, 0, logUnit)
    call UTCounter_startCounter(counterTarget)
    isTargetMet = .false.
    do while (.not.isTargetMet)
       call UTCounter_incrementCounter(2)
       call UTCounter_progressCounter(isTargetMet)
    end do

    call UTCounter_stopCounter(doAsyncReturn=.true.)
    call UTCounter_finalize()

    close(logUnit)
  end subroutine test3


  !Test an early shut-down before the counter target is met
  subroutine test4()
    implicit none
    integer, parameter :: logUnit = 24
    character(len=100) :: fileName
    integer :: counterTarget, errorcode, ierr, i
    logical :: isTargetMet

    counterTarget = mpisize + 1 !One higher than the number of MPI ranks
    if (logUnit /= 6) then
       write(fileName,'(a,i6.6,a)') 'test4_', mpirank, '.out'
       open (logUnit, file=fileName)
    end if

    call UTCounter_init(mpiComm, mpiMasterRank, logUnit)
    call UTCounter_startCounter(counterTarget)

    call UTCounter_incrementCounter(1)
    do i = 1, 1000
       call UTCounter_progressCounter(isTargetMet)
    end do

    if (isTargetMet) then
       print *, "MPI rank:", mpirank, &
            ", ERROR!!!  It should not be possible to reach the target"
       if (logUnit /= 6) then
          close(logUnit)
       end if
       call MPI_Abort(mpiComm, errorcode, ierr)
       stop
    end if

    !UTCounter_stopCounter forces completion/cancelation of all
    !outstanding commnunications.  There are still unmatched
    !communications because the counter target was not met.  This
    !means we must synchronize so that the outstanding communications
    !do not interfere with the next unit test.
    call UTCounter_stopCounter(doAsyncReturn=.false.)
    call UTCounter_finalize()

    if (logUnit /= 6) then
       close(logUnit)
    end if
  end subroutine test4
  


  !Everyone contributes 1 to the total count
  subroutine test5()
    implicit none
    integer, parameter :: logUnit = 25
    integer :: counterTarget
    character(len=100) :: fileName
    logical :: isTargetMet

    counterTarget = mpiSize
    write(fileName,'(a,i6.6,a)') 'test5_', mpirank, '.out'
    open (logUnit, file=fileName)

    call UTCounter_init(mpiComm, mpiMasterRank, logUnit)
    call UTCounter_startCounter(counterTarget)

    call UTCounter_incrementCounter(1)
    isTargetMet = .false.
    do while (.not.isTargetMet)
       call UTCounter_progressCounter(isTargetMet)
    end do

    !No need to synchronize because we can guarantee that
    ! 1. Communication is done
    ! 2. There are no outstanding contribution messages from the
    ! slave MPI ranks.
    call UTCounter_stopCounter(doAsyncReturn=.true.)
    call UTCounter_finalize()
    
    close(logUnit)
  end subroutine test5

 
  !Ranks contribute an unspecified count
  subroutine test6()
    implicit none
    integer, parameter :: logUnit = 26
    integer :: counterTarget, finalCount
    logical :: isTargetMet
    character(len=100) :: fileName

    counterTarget = 500*(mpiSize + 1)
    write(fileName,'(a,i6.6,a)') 'test6_', mpirank, '.out'
    open (logUnit, file=fileName)

    call UTCounter_init(mpiComm, mpiMasterRank, logUnit)
    call UTCounter_startCounter(counterTarget)

    isTargetMet = .false.
    do while (.not.isTargetMet)
       call UTCounter_incrementCounter(1)
       call UTCounter_progressCounter(isTargetMet, finalCount)
    end do
    write(logUnit,*) "Target was ", counterTarget, &
         ". Final count is ", finalCount

    !Must finalize in the standard way, i.e. with synchronization,
    !because there are possibly outstanding slave contributions
    call UTCounter_stopCounter()
    call UTCounter_finalize()

    close(logUnit)
  end subroutine test6



  !This is THE stress test which emulates the usage of the UTCounter
  !module in the laser.  The cumulative sum of all contribution
  !messages adds up *exactly* to the counter target.  No
  !synchronization is needed because there are no outstanding sends
  !that can be intercepted in different invocations of the UTCounter
  !module.
  subroutine test7()
    implicit none
    integer, parameter :: logUnit = 27, countsPerRank = 5, &
         counterInvocations = 100
    character(len=100) :: fileName
    integer :: counterTarget, i
    logical :: isTargetMet

    counterTarget = mpisize * countsPerRank
    write(fileName,'(a,i6.6,a)') 'test7_', mpirank, '.out'
    open (logUnit, file=fileName)

    call UTCounter_init(mpiComm, mpiMasterRank, logUnit)

    do i = 1, counterInvocations
       call UTCounter_startCounter(counterTarget)
       call UTCounter_incrementCounter(countsPerRank)

       isTargetMet = .false.
       do while (.not.isTargetMet) !Spin until everyone is done
          call UTCounter_progressCounter(isTargetMet)
       end do

       !No need to synchronize because we can guarantee that
       ! 1. Communication is done
       ! 2. There are no outstanding contribution messages from the
       ! slave MPI ranks.
    end do
    call UTCounter_finalize()

    close(logUnit)
  end subroutine test7

end subroutine test_subroutine
