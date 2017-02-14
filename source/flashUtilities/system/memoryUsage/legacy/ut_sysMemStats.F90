!!****if* source/flashUtilities/system/memoryUsage/legacy/ut_sysMemStats
!!
!! NAME
!!  ut_sysMemStats
!!
!! SYNOPSIS
!!
!!  ut_sysMemStats(integer(IN)             :: verbosity,
!!                 integer(IN)             :: memorySamplerArg,
!!                 type(meminfo_t)(:)(OUT) :: meminfo
!!                 integer(OUT)            :: numStats)
!!  
!! DESCRIPTION
!!  
!!  Returns memory usage information.
!!
!! ARGUMENTS
!!
!!  verbosity : The verbosity of the memory usage information.  There are
!!              currently two levels
!!               0 = regular
!!               1 = more verbose
!!
!!  memorySamplerArg : The method for sampling memory.  The following are valid
!!                   UT_SYSMEM_PROC = use /proc/self/stat
!!                   UT_SYSMEM_RUSAGE = use rusage
!!                   UT_SYSMEM_MALLINFO = use mallinfo (f2003)
!!                   UT_SYSMEM_BGKERNEL = use Kernel_GetMemorySize (f2003)
!!                   UT_SYSMEM_AUTO = use the default sampler for this machine
!!                   UT_SYSMEM_ALL = use all samplers supported by this machine
!!                  (f2003) indicates that the sampler is only supported
!!                  by the memory usage code in directory f2003.  This is
!!                  only included in a FLASH application if you setup with
!!                  useFortran2003=True.
!!
!!  meminfo : An array containing the returned memory information.
!!
!!  numStats : The number of memory statistics in meminfo.
!!
!!***

#include "constants.h"
#include "ut_sysMem.h"

subroutine ut_sysMemStats(verbosity, memorySamplerArg, meminfo, numStats)
  use ut_sysMemData, ONLY : meminfo_t
  implicit none
  integer, intent(IN) :: verbosity, memorySamplerArg
  type(meminfo_t), dimension(:), intent(OUT) :: meminfo
  integer, intent(OUT) :: numStats

  real :: vsizeMB, resMemMB
  integer :: meminfoSize, ierr, memSampler
  integer, external :: ut_sys_mem_usage

  meminfoSize = size(meminfo,1)
  ierr = ut_sys_mem_usage(vsizeMB, resMemMB, memSampler)
  numStats = 0  

  select case (memSampler)
  case (UT_SYSMEM_PROC)
     if (meminfoSize >= 2) then
        numStats = 2
        meminfo(1) % description = "/proc vsize    (MB):"
        meminfo(1) % measurement = vsizeMB
        meminfo(2) % description = "/proc rss      (MB):"
        meminfo(2) % measurement = resMemMB
     end if
  case (UT_SYSMEM_RUSAGE)
     if (meminfoSize >= 1) then
        numStats = 1
        meminfo(1) % description = "rusage rss     (MB):"
        meminfo(1) % measurement = resMemMB
     end if
  end select
end subroutine ut_sysMemStats
