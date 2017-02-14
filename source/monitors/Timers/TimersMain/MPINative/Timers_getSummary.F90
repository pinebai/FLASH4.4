!!****if* source/monitors/Timers/TimersMain/MPINative/Timers_getSummary
!!
!! NAME
!!   
!!  Timers_getSummary
!!
!! SYNOPSIS
!!
!!
!!  Timers_getSummary(integer(in) :: nIntervals)
!!
!! DESCRIPTION
!!
!!  Write out the Timers summary to the logfile.  Timers collects the
!!  performance data and then calls Logfile_writeSummary to do the
!!  actual formatting.
!!
!!  It should be safe to call this at any point in the simulation-- an
!!  updated summary will appear in the logfile.
!!
!!  The summary written only includes data from the master processor.
!!
!! ARGUMENTS
!!
!!  nIntervals - number of subintervals timed, which is determined 
!!               by the caller of this code, but this will 
!!               typically be the number of timesteps taken
!!               since the last time Timers_init was called. 
!!
!! NOTES 
!!
!!  The ability to count the number of zones evolved is
!!  currently not supported in this implementation
!!
!!***

! Requesting diagnostic output when timers on different procs are inconsistent:
#define DEBUG_TIMERSDIFFS
! Uncomment to get more diagnostic output (at end of each run):
!#define DEBUG_TIMERS

subroutine Timers_getSummary(nIntervals)

  use Timers_data, ONLY: tmr_stack, tmr_bigInt, tmr_acctSegs, &
       tmr_numsegments, tmr_maxtimerparents, tmr_initDate, tmr_initTime, &
       tmr_writeStatSummary, tmr_eachProcWritesSummary, tmr_globalMe

  use Logfile_interface, ONLY : Logfile_writeSummary, Logfile_stamp
  use Grid_interface, ONLY : Grid_getLocalNumBlks, &
    Grid_getBlkIndexLimits

  implicit none

#include "constants.h"

  integer, intent(in) ::  nIntervals
  integer, parameter :: maxRows = 1000
  integer, parameter :: maxColumns = 10
  character(len=MAX_STRING_LENGTH), dimension(maxRows, maxColumns) :: perfmonArr
  
  integer             :: length, dim
  integer             :: index
  integer             :: totalSimSteps, globalNumBlocks
  integer             :: i, j
  real                :: totalWCTime !wall clock time
  real                :: zonesPerSecond
  real                :: curTime
  integer             :: numHeaders
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  character(len=40), save      :: dateStr
  character(len=MAX_STRING_LENGTH) :: numToStr
  integer (kind=tmr_bigInt)           :: numEvolvedZones
  type(tmr_stack) :: rootStack

  logical :: allSame

  real :: tempTime

  !This initialization is here so that TAU generates valid code.
  !Without it, TAU declares a variable after !$omp master.
  tempTime = 0.0

#ifdef NOOP
  return
#endif

  !$omp master
  call tmr_stackZero(rootStack)

  totalSimSteps = nIntervals

  call tmr_etime(curTime)
  totalWCTime = curTime - tmr_initTime
  call current_date_time(dateStr)
  
!!$  !old DEV: code for getting an estimate of 'zones per second' measure - disabled.
!!$  !old DEV: this would have to become more complicated with non fixed block size
!!$  call Grid_getLocalNumBlks(globalNumBlocks)
!!$  call Grid_getBlkIndexLimits(1, blkLimits, blkLimitsGC)
!!$  numEvolvedZones = globalNumBlocks
!!$  do i=1, MDIM
!!$     numEvolvedZones = numEvolvedZones * (blkLimits(HIGH,i) - blkLimits(LOW,i) + 1)
!!$  end do
!!$  
!!$  zonesPerSecond = numEvolvedZones / totalWCTime
  
  do i = 1, tmr_numSegments
     do j = 1, tmr_maxTimerParents
        if (tmr_acctSegs(i)%isTimed(j)) then
           call tmr_etime(tempTime)
           tmr_acctSegs(i)%time(j)       = tempTime - tmr_acctSegs(i)%dtime(j) + tmr_acctSegs(i)%time(j)
           tmr_acctSegs(i)%dtime(j)      = 0.
           tmr_acctSegs(i)%isTimed(j) = .false.
        endif
        tmr_acctSegs(i)%pctTime(j) = 100. * tmr_acctSegs(i)%time(j) / totalWCTime
     enddo
  enddo
  
  write (perfmonArr(1,1), "(A)") 'beginning'
  write (perfmonArr(1,2), "(A)") trim(tmr_initDate)
  
  write (perfmonArr(2,1), "(A)") 'ending'
  write (perfmonArr(2,2), "(A)") trim(dateStr)
  
  write (perfmonArr(3,1), "(A)") 'seconds in monitoring period'
  write (perfmonArr(3,2), "(F15.3)") totalWCTime
  
  write (perfmonArr(4,1), "(A)") 'number of evolution steps'
  write (perfmonArr(4,2),     *) totalSimSteps
  
!!  these numbers are not accurate at all, so we don't include 
!!  them in output

!!$  !code for output of estimate of 'zones per second' measure - disabled.
!!$  write (perfmonArr(5,1), "(A)") 'number of evolved zones'
!!$  write (numToStr, "(I15)") numEvolvedZones
!!$  write (perfmonArr(5,2),  "(A)") trim(adjustl(numToStr))
!!$
!!$  write (perfmonArr(6,1), "(A)") 'zones per second'
!!$  write (perfmonArr(6,2), "(F20.3)") zonesPerSecond
  
! ==== Prepare original summary that shows data from the master proc only ====
  write (perfmonArr(5,2), "(A)") 'accounting unit'
  write (perfmonArr(5,3), "(A)") 'time secs'
  write (perfmonArr(5,4), "(A)") 'num calls'
  write (perfmonArr(5,5), "(A)") 'secs avg'
  write (perfmonArr(5,6), "(A)") 'time pct'
  
  numHeaders = 4
  index = 6
  call tmr_buildSummary(perfmonArr, maxRows, maxColumns, index, 0, rootStack, .FALSE.)
  
  length = index-1
  dim = 6
  
  ! always write process zero summary to the logfile
  call Logfile_writeSummary(perfmonArr(1:length,1:dim), length, dim, MAX_STRING_LENGTH, numHeaders)        
  ! now each process write summary to its own file
  if (tmr_eachProcWritesSummary) then
     call Logfile_writeSummary(perfmonArr(1:length,1:dim), length, dim, MAX_STRING_LENGTH, numHeaders, separateFiles=.true.)        
  end if
  
! ==== Prepare additional summary that does combine data from several processors ====
! ==== But decide to write it based on a runtime parameter and a check to make sure 
! ==== the write won't hang because of different timer stacks ===

  ! print a warning if rt parm tells us not to write stats
  if (.not.tmr_writeStatSummary) then
     call Logfile_stamp( &
          'Not writing timer max/min/avg values because runtime parameter writeStatSummary == .false.', &
          '[Timers_getSummary]')
  else
     
     ! code will hang if not all processors have the same timers; also,
     ! doing the summary if they're not the same would produce garbage

     allSame = .true.
     
     call tmr_allProcsSameTimers(tmr_globalMe, allSame)

     if (.not.allSame) then
        call Logfile_stamp( &
             'Not writing timer max/min/avg values because not all processors had same timers', &
             '[Timers_getSummary]')
     else
        write (perfmonArr(5,3), "(A)") 'max/proc (s)'
        write (perfmonArr(5,4), "(A)") 'min/proc (s)'
        write (perfmonArr(5,5), "(A)") 'avg/proc (s)'
        write (perfmonArr(5,6), "(A)") 'num calls'
        
        numHeaders = 4
        index = 6
        call tmr_buildSummary(perfmonArr, maxRows, maxColumns, index, 0, rootStack, .TRUE.)
  
        length = index-1
        dim = 6
     
        ! now write the reduced summary to the logfile
        call Logfile_writeSummary( perfmonArr(1:length,1:dim), length, dim, &
             MAX_STRING_LENGTH, numHeaders,reduced=.TRUE.,separateFiles=.FALSE.)
     end if
  end if
  !$omp end master
  return

end subroutine Timers_getSummary


subroutine tmr_allProcsSameTimers(myPE, areSame)
  
  use Timers_data, ONLY: tmr_acctSegType, tmr_acctSegs, tmr_maxSegments, &
       tmr_numSegments, tmr_globalComm, tmr_globalMe

  implicit none

#include "Flash_mpi.h"

  integer, intent(in)  :: myPE
  logical, intent(out) :: areSame


  type (tmr_acctSegType) :: rootSegs(tmr_maxSegments)

  integer :: i, j
  logical :: stackListsSame
  integer :: numRootTimers
  logical :: numTimersAllSame, areSameReduction
  integer :: ierr

  areSame = .true.
  
  ! first, see if number of timers is the same
  numRootTimers = tmr_numSegments
  call mpi_bcast(numRootTimers, 1, FLASH_INTEGER, MASTER_PE, tmr_globalComm, ierr)
  call mpi_allreduce(numRootTimers .eq. tmr_numSegments, numTimersAllSame, 1, MPI_LOGICAL, &
       MPI_LAND, tmr_globalComm, ierr)

  if (.not.numTimersAllSame) then
     areSame = .false.
#ifdef DEBUG_TIMERSDIFFS

#ifdef DEBUG_TIMERDIFFS_DETAILED
     print*,'ERROR NUM TIMERS: ', tmr_globalMe, tmr_numSegments

     if(tmr_globalMe == MASTER_PE .or. numRootTimers /= tmr_numSegments) then
        do i = 1, tmr_numSegments
           print '(i10,i6,2x,a)', tmr_globalMe, i, trim(tmr_acctSegs(i)%name)
        end do
     end if

#else
     print*,'proc',tmr_globalMe,' says (.not.numTimersAllSame), whaddayamakeofthat?'  !DEBUG
#endif

#endif
  else
     ! communicate root timers to everyone
     call tmr_broadcastRootTimers(tmr_globalMe, numRootTimers, rootSegs)

     if (tmr_globalMe .NE. MASTER_PE) then
     ! compare names  and stacks of my timers to the ones I received from root;
     ! beware, if the data is not correct in rootSegs after the transfer,
     ! the following stackList comparison operations might go out of bounds
        do i = 1, numRootTimers
           if (rootSegs(i)%name .ne. tmr_acctSegs(i)%name) then
              areSame = .false.
#ifdef DEBUG_TIMERSDIFFS
              print*,'proc',tmr_globalMe,' says about root timer',i,': (',rootSegs(i)%name,'.ne.',tmr_acctSegs(i)%name,').'  !DEBUG

#ifdef DEBUG_TIMERDIFFS_DETAILED
              do j = 1, numRootTimers
                 print '(i10,i5,2x,a30,2x,a30)', tmr_globalMe, j, trim(rootSegs(j)%name), trim(tmr_acctSegs(j)%name)
              end do
#endif
#endif
              exit
           else
              ! now make sure my call stacks are the same as root's
              call tmr_stackListsEqual(tmr_acctSegs(i)%stacks, rootSegs(i)%stacks, stackListsSame)
              if (.not. stackListsSame) then
                 areSame = .false.
#ifdef DEBUG_TIMERSDIFFS
                 print*,'proc',tmr_globalMe,' says about timer',i,': (.not.stackListsSame).'  !DEBUG
#endif
                 exit
              end if
           end if
        end do
     end if
  end if

  call mpi_allreduce(areSame, areSameReduction, 1, MPI_LOGICAL, &
       MPI_LAND, tmr_globalComm, ierr)

  areSame = areSameReduction
  
  return
end subroutine tmr_allProcsSameTimers
    
subroutine tmr_broadcastRootTimers(mype, numRootTimers, rootSegs)
  
  use Timers_data, ONLY: tmr_acctSegType, tmr_acctSegs, tmr_nameSize, &
       tmr_nameSize, tmr_maxTimerParents, tmr_maxCallStackDepth, &
       tmr_maxSegments,tmr_globalComm, tmr_globalMe

  implicit none

#include "Flash_mpi.h"

  integer, intent(in) :: mype
  integer, intent(in) :: numRootTimers
  type(tmr_acctSegType), intent(out) :: rootSegs(tmr_maxSegments)

  integer :: numPEs
  integer :: ierr
  integer :: charExtent, realExtent, logicalExtent, integerExtent
  integer :: acctSegExtent
  integer :: count, blockcounts(4), offsets(5), oldtypes(4), acctSegMpiType
  integer :: numberOfIntegersInAcctSeg
  integer :: i

  integer :: status(MPI_STATUS_SIZE), iSendReq
  integer, allocatable :: statuses(:,:)
  integer, allocatable :: request(:)

  call MPI_COMM_SIZE(tmr_globalComm, numPEs, ierr)

  allocate(request(numPEs))

  call MPI_Type_extent(FLASH_REAL, realExtent, ierr)
  call MPI_Type_extent(FLASH_INTEGER, integerExtent, ierr)
  call MPI_Type_extent(MPI_CHARACTER, charExtent, ierr)
  call MPI_Type_extent(MPI_LOGICAL, logicalExtent, ierr)

  ! setup the description of the mpi type for an accounting segment

  ! setup the names field in the mpi type
  offsets(1) = 0
  call MPI_ADDRESS(tmr_acctSegs(1)%name, offsets(1), ierr)
  oldtypes(1) = MPI_CHARACTER
  blockcounts(1) = tmr_nameSize

  ! setup the reals
  offsets(2) = tmr_nameSize * charExtent
  call MPI_ADDRESS(tmr_acctSegs(1)%time, offsets(2), ierr)
  oldtypes(2) = FLASH_REAL
  blockcounts(2) = tmr_maxTimerParents * 4

  ! setup the logicals
  offsets(3) = offsets(2) + tmr_maxTimerParents * 4 * realExtent
  call MPI_ADDRESS(tmr_acctSegs(1)%isTimed, offsets(3), ierr)
  oldtypes(3) = MPI_LOGICAL
  blockcounts(3) = tmr_maxTimerParents
  
  ! setup the integers
  offsets(4) = offsets(3) + tmr_maxTimerParents * logicalExtent
  call MPI_ADDRESS(tmr_acctSegs(1)%timesCalled, offsets(4), ierr)
  oldtypes(4) = FLASH_INTEGER
  numberOfIntegersInAcctSeg = tmr_maxTimerParents + &  ! this includes timesCalled and stacks fields
       2 + tmr_maxTimerParents * (tmr_maxCallStackDepth + 1) 
  blockcounts(4) = numberOfIntegersInAcctSeg


  call MPI_ADDRESS(tmr_acctSegs(2)%name, offsets(5), ierr)
  offsets(2) = offsets(2) - offsets(1)
  offsets(3) = offsets(3) - offsets(1)
  offsets(4) = offsets(4) - offsets(1)
  offsets(5) = offsets(5) - offsets(1)
  offsets(1) = 0

  call MPI_TYPE_STRUCT(4, blockcounts, offsets, oldtypes, &
       acctSegMpiType, ierr)
  call MPI_TYPE_COMMIT(acctSegMpiType, ierr)
  call MPI_TYPE_EXTENT(acctSegMpiType, acctSegExtent, ierr)

#ifdef DEBUG_TIMERS
  if (tmr_globalMe .eq. MASTER_PE) then
     print*,'tmr_broadcastRootTimers: numRootTimers=',numRootTimers  !DEBUG
     print*,'max. lastItem:',maxval(tmr_acctSegs%stacks%lastItem) !DEBUG
!     print*,'sizeof(tmr_acctSegs(1)),MPI_TYPE_EXTENT(acctSegMpiType),SUM(extents),offset(tmr_acctSegs(2)) ='
!     print*, sizeof(tmr_acctSegs(1)), acctSegExtent,&
     print*,'MPI_TYPE_EXTENT(acctSegMpiType),SUM(extents),offset(tmr_acctSegs(2)) ='
     print*,  acctSegExtent,&
          tmr_nameSize * charExtent  +  tmr_maxTimerParents * 4 * realExtent  +  tmr_maxTimerParents * logicalExtent  + &
          numberOfIntegersInAcctSeg * integerExtent, &
          offsets(5)
  end if
#endif

#define USE_MPI_BCAST_HERE

#ifdef USE_MPI_BCAST_HERE
  rootSegs(1:numRootTimers) = tmr_acctSegs(1:numRootTimers)
!  call MPI_SENDRECV(tmr_acctSegs, numRootTimers, acctSegMpiType, MASTER_PE,0, &
!       rootSegs, numRootTimers, acctSegMpiType, MASTER_PE,0, tmr_globalComm, ierr)
  call MPI_BCAST(rootSegs(1)%name, numRootTimers, acctSegMpiType, MASTER_PE, tmr_globalComm, ierr)
#else
  if (tmr_globalMe .ne. MASTER_PE) then
     call MPI_IRECV(rootSegs, numRootTimers, acctSegMpiType, MASTER_PE, tmr_globalMe, tmr_globalComm, request(tmr_globalMe+1), ierr)
     if (ierr .ne. MPI_SUCCESS) print*,'ierr is',ierr,' after MPI_IRECV in tmr_broadcastRootTimers!'
  end if

  ! doing this with a regular mpi_send can deadlock if there are enough timers because enough timers
  ! make the message big enough to not fit in the standard buffer, so it waits for the posted receive...
  ! root posts its receive after its send, hence the hang
  if (tmr_globalMe .eq. MASTER_PE) then
     iSendReq = 0
     do i = 0, numPEs-1
        if (i .ne. MASTER_PE) then
           iSendReq = iSendReq + 1
           call MPI_ISEND(tmr_acctSegs, numRootTimers, acctSegMpiType, i, i, tmr_globalComm, request(iSendReq), ierr)
           if (ierr .ne. MPI_SUCCESS) print*,'ierr is',ierr,' after MPI_ISEND in tmr_broadcastRootTimers!'
        end if
     end do
  end if

  if (tmr_globalMe .ne. MASTER_PE) then
     call MPI_WAIT(request(tmr_globalMe+1), status, ierr)
     if (ierr .ne. MPI_SUCCESS) print*,'ierr is',ierr,' after MPI_WAIT in tmr_broadcastRootTimers!'
  else
     allocate(statuses(MPI_STATUS_SIZE,iSendReq))
     call MPI_WAITALL(isendreq,request, statuses, ierr)
     if (ierr .ne. MPI_SUCCESS) print*,'ierr is',ierr,' after MPI_WAITALL in tmr_broadcastRootTimers!'
     deallocate(statuses)
  end if
#endif

  call MPI_TYPE_FREE(acctSegMpiType, ierr)
  deallocate(request)

end subroutine tmr_broadcastRootTimers
