Program test_pipeline
  call test_subroutine
End Program test_pipeline


subroutine test_subroutine
  use mpi
  use UTPipeline
  implicit none
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
  if (mpiRank == 0) print *, "starting test2"
  call test2()
  if (mpiRank == 0) print *, "finished all tests!!"

  call MPI_Errhandler_set(mpiComm, MPI_ERRORS_ARE_FATAL, ierr)
  call MPI_Finalize(ierr)

contains

  !Test a simple start-up and shut-down
  subroutine test1()
    implicit none
    integer :: itemSize = 1, maxItems = 1, channelSize = 1, numChannels = 0
    integer, dimension(1) :: procList
    call UTPipeline_init(itemSize, maxItems, channelSize, mpiComm, numChannels, procList)
    call UTPipeline_finalize()
  end subroutine test1


  subroutine test2()
    implicit none
    integer :: mpiCartComm, ierr, mpiCartRank, mpiCartSize, &
         lneigh, rneigh, numChannels, numItems, channelSize, i
    integer, parameter :: ndims = 1, direction = 0, disp = 1
    logical, parameter :: reorder = .false.
    integer, dimension(ndims) :: dim_size
    logical, dimension(ndims) :: periods
    integer, allocatable, dimension(:) :: procList
    integer, parameter :: itemSize = 1, maxItems = 1, hop_target = 10
    real, dimension(itemSize,maxItems) :: items
    logical :: isHandled, isClosing, isDone, doFlush, localDone
    character(len=100) :: fileName
    integer, parameter :: logUnit = 24
    dim_size(1) = mpiSize
    periods(1) = .true.


    channelSize = maxItems !The channel Size must be >= maxItems in this test.


    call MPI_Cart_create(mpiComm, ndims, dim_size, periods, reorder, &
         mpiCartComm, ierr)
    call MPI_Comm_size(mpiCartComm, mpiCartSize, ierr)
    call MPI_Comm_rank(mpiCartComm, mpiCartRank, ierr)
    call MPI_Cart_shift(mpiCartComm, direction, disp, lneigh, rneigh, ierr)
    print *, "Rank", mpiCartRank, "lneigh", lneigh, "rneigh", rneigh

    if (lneigh == rneigh) then
       allocate(procList(1))
       procList(1) = lneigh
       numChannels = 1
    else
       allocate(procList(2))
       procList(1) = lneigh
       procList(2) = rneigh
       numChannels = 2
    end if

    write(fileName,'(a,i6.6,a)') 'test2_', mpiCartRank, '.out'
    open (logUnit, file=fileName)


    items(:,:) = 0
    call UTPipeline_init(itemSize, maxItems, channelSize, mpiCartComm, numChannels, procList, logUnit)
    call UTPipeline_initComm()


    numItems = maxItems
    doFlush = .true. !Needed when channelSize > maxItems
    localDone = .false.
    do while (.not.localDone)
       !Everyone send to the right.  Fill up the internal pipeline buffers
       do i = 1, numItems
          call UTPipeline_sendItem(items(:,1), rneigh, isHandled)
          if (.not.isHandled) then
             call Driver_abortFlash('Test arranged so sends always handled')
          end if
       end do

       !Everyone receive
       numItems = 0
       do while (numItems < maxItems)
          call UTPipeline_progressComm(doFlush)
          call UTPipeline_getItems(items, maxItems, numItems)
       end do
       
       !Update the hop count
       items(1,1:numItems) = items(1,1:numItems) + 1
       localDone = all(items(1,1:numItems) > hop_target)
    end do

    isClosing = .false.
    do while (.not.isClosing)
       call UTPipeline_closeSendChannels(isClosing)
       call UTPipeline_progressComm()
    end do

    isDone = .false.
    do while (.not.isDone)
       call UTPipeline_isDone(isDone)
       call UTPipeline_progressComm()
    end do

    call UTPipeline_finalizeComm()
    call UTPipeline_finalize()

    deallocate(procList)

    close(logUnit)

    call MPI_Comm_free(mpiCartComm, ierr)
  end subroutine test2

end subroutine test_subroutine
