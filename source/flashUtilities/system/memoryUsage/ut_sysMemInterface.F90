!!****ih* source/flashUtilities/system/memoryUsage/ut_sysMemInterface
!!
!! NAME
!!  ut_sysMemInterface 
!!
!! SYNOPSIS
!!  use ut_sysMemInterface
!!
!!***

#include "constants.h"

module ut_sysMemInterface
  implicit none
  interface
     subroutine ut_sysMemSummaryStats(comm, verbosity, memorySampler, &
          memsummary, numStats)
       use ut_sysMemData, ONLY : memsummary_t, meminfo_t
       implicit none
       integer, intent(IN) :: comm, verbosity, memorySampler
       type(memsummary_t), dimension(:), intent(OUT) :: memsummary
       integer, intent(OUT) :: numStats
     end subroutine ut_sysMemSummaryStats
  end interface

  interface
     subroutine ut_sysMemStats(verbosity, memorySamplerArg, meminfo, numStats)
       use ut_sysMemData, ONLY : meminfo_t
       implicit none
       integer, intent(IN) :: verbosity, memorySamplerArg
       type(meminfo_t), dimension(:), intent(OUT) :: meminfo
       integer, intent(OUT) :: numStats
     end subroutine ut_sysMemStats
  end interface
end module ut_sysMemInterface
