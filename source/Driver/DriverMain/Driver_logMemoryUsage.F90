!!****if* source/Driver/DriverMain/Driver_logMemoryUsage
!!
!! NAME
!!
!!  Driver_logMemoryUsage
!!
!! SYNOPSIS
!!
!!  Driver_logMemoryUsage(character(len=*)(IN) :: callsite)
!!
!! DESCRIPTION
!!
!!  Logs memory usage
!!
!!
!! ARGUMENTS
!!
!!  callsite -    A string storing from where we call this subroutine.
!!
!! NOTES
!!
!!  This routine is a no-op if the preprocessor symbol FLASH_USE_MEMORYUSAGE is undefined.
!!***


subroutine Driver_logMemoryUsage (callsite)

#include "constants.h"
#include "Flash.h"

#ifdef FLASH_USE_MEMORYUSAGE
#include "ut_sysMem.h"
  use Driver_data, ONLY : dr_globalMe, dr_globalComm
  use ut_sysMemInterface, ONLY : ut_sysMemSummaryStats
  use ut_sysMemData, ONLY : memsummary_t
  use Logfile_interface, ONLY: Logfile_stamp
#endif
  implicit none
  character(len=*), intent(IN) :: callsite
#ifdef FLASH_USE_MEMORYUSAGE
  integer, parameter :: maxStats = 20, verbosity = 1, &
       memorySampler = UT_SYSMEM_ALL
  type(memsummary_t), dimension(maxStats) :: memSummary
  integer :: i, numStats, sLen
  character (len=12) :: minMem_string, maxMem_string, avgMem_string
  character (len=12) :: minMemRank_string, maxMemRank_string
  character (len=MAX_STRING_LENGTH) :: rankStr

  call ut_sysMemSummaryStats(dr_globalComm, verbosity, memorySampler, &
       memSummary, numStats)

  if (dr_globalMe == MASTER_PE) then
     call Logfile_stamp(callsite, tag="memory call site")
     do i = 1, numStats
        write (minMem_string, '(f12.2)') memSummary(i) % min
        write (maxMem_string, '(f12.2)') memSummary(i) % max
        write (avgMem_string, '(f12.2)') memSummary(i) % avg
        write (minMemRank_string, '(i12)') memSummary(i) % min_rank
        write (maxMemRank_string, '(i12)') memSummary(i) % max_rank

        call Logfile_stamp(trim(memSummary(i) % description) // &
             trim(minMem_string) // ' (min)  ' // &
             trim(maxMem_string) // ' (max)  ' // &
             trim(avgMem_string) // ' (avg) ', tag="memory")

        if (verbosity > 0) then
           sLen = len_trim(memSummary(i) % description)
           rankStr = ''
           if (sLen > 0) then
              rankStr = repeat(' ',sLen)
              if (sLen > 4) then
                 rankStr(sLen-4:sLen) = 'rank:'
              else
                 rankStr(sLen:sLen) = ':'
              end if
           end if

           call Logfile_stamp(trim(rankStr) // &
                trim(minMemRank_string) // ' (min)  ' // &
                trim(maxMemRank_string) // ' (max)  ', tag="memory")
        end if
     end do
  end if
#endif
end subroutine Driver_logMemoryUsage
