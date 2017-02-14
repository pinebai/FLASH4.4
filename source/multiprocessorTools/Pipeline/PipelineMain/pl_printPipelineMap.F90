!!****if* source/multiprocessorTools/Pipeline/PipelineMain/pl_printPipelineMap
!!
!! NAME
!!
!!  pl_printPipelineMap
!!
!! SYNOPSIS
!!
!!  call pl_printPipelineMap ()
!!
!! DESCRIPTION
!!
!!  Prints the pipeline map to the master's log file (if log files are requested).
!!
!! ARGUMENTS
!!
!!  none
!!
!! NOTES
!!
!!  Only the master rank can do the printing as he is the only one who has the
!!  pipeline map.
!!
!!***

subroutine pl_printPipelineMap ()

  use Pipeline_data, ONLY : pl_doLog,       &
                            pl_logUnit,     &
                            pl_masterRank,  &
                            pl_maxChannels, &
                            pl_pipelineMap, &
                            pl_rank,        &
                            pl_size

  implicit none

  character (len=6),  parameter :: empty6    = '      '
  character (len=6),  parameter :: line6     = '------'
  character (len=18), parameter :: line18    = '-----------------|'
  character (len=18), parameter :: rankLabel = '         Rank -> |'
  character (len=18), parameter :: row1Label = '   # of Channels |'
  character (len=18), parameter :: row2Label = ' Channel proc ID |'

  character (len=6) :: lineOfMap (1:10)

  logical :: printDone

  integer :: channel
  integer :: n, number
  integer :: rank
  integer :: rankBeg, rankEnd, rankMax, rankCnt

  integer, parameter :: nRanksPerBlock = 10
!
!
!     ...Immediate return if rank is not master rank or no logfiles are
!        requested.
!
!
  if (.not.pl_doLog .or. pl_rank /= pl_masterRank) then
       return
  end if
!
!
!     ...Print out the title.
!
!
  write (pl_logUnit,'(/)')
  write (pl_logUnit,'(a)') ' Pipeline Map '
  write (pl_logUnit,'(/)')
!
!
!     ...Print the pipeline map in blocks of ranks.
!
!
  rankBeg = 0
  rankMax = pl_size - 1
  printDone = .false.

  do while (.not.printDone)

     rankEnd = min (rankMax, rankBeg + nRanksPerBlock - 1)
     rankCnt = rankEnd - rankBeg + 1

     write (pl_logUnit,'(a,10i6)') rankLabel,  (rank , rank  =  rankBeg,rankEnd)
     write (pl_logUnit,'(11(a))') line18, (line6, rank = 1,rankCnt) 
     write (pl_logUnit,'(a,10i6)') row1Label, pl_pipelineMap (1,rankBeg:rankEnd)
     write (pl_logUnit,'(11(a))') line18, (line6, rank = 1,rankCnt) 

     do channel = 1, pl_maxChannels

        n = 0
        do rank = rankBeg, rankEnd
           n = n + 1
           number = pl_pipelineMap (1+channel,rank)
           if (number /= -1) then
               write (lineOfMap (n),'(i6)') number
           else
               lineOfMap (n) = empty6     ! do not print non-existing channels with label -1
           end if
        end do

        write (pl_logUnit,'(a,10a6)') row2Label, lineOfMap (1:n)
     end do

     write (pl_logUnit,'(/)')

     rankBeg = rankEnd + 1
     printDone = rankBeg > rankMax

  end do
!
!
!     ...Ready!
!
!
  return
end subroutine pl_printPipelineMap
