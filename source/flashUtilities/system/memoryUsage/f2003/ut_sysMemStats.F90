!!****if* source/flashUtilities/system/memoryUsage/f2003/ut_sysMemStats
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
  use iso_c_binding, only : C_NULL_CHAR, c_ptr, c_char, c_f_pointer, c_associated, &
       c_f_procpointer, c_funptr, c_funloc, C_NULL_PTR, c_loc, C_NULL_FUNPTR
  use ut_sysMemCInterface, ONLY : ut_sysMem, ut_sysMemNumStats, ut_sysMemMallinfo, &
       ut_sysMemMallinfoNumStats, ut_sysMemProc, ut_sysMemProcNumStats, &
       ut_sysMemBGKernel, ut_sysMemBGKernelNumStats, &
       ut_sysMemRusage, ut_sysMemRusageNumStats, ut_sysMemAutoDetect
  use ut_sysMemData, ONLY : meminfo_t
  use ut_sysMemCData, ONLY : c_meminfo_t
  implicit none
  integer, intent(IN) :: verbosity, memorySamplerArg
  type(meminfo_t), dimension(:), intent(OUT) :: meminfo
  integer, intent(OUT) :: numStats

  type(c_meminfo_t), dimension(size(meminfo,1)) :: c_meminfo
  integer, parameter :: AllFns = 4
  type(c_funptr), dimension(AllFns) :: c_num_stats_fn, c_memory_usage_fn
  procedure(ut_sysMemNumStats), pointer :: num_stats_fn
  procedure(ut_sysMem), pointer :: memory_usage_fn
  integer :: i, j, m, freeSlot, fnStats, indexPtr, numFns, meminfoSize, memorySampler
  character(kind=c_char), dimension(:), pointer :: msgPtr
  character(len=MAX_STRING_LENGTH) :: description
  character :: c
  meminfoSize = size(meminfo,1)

  if (memorySamplerArg == UT_SYSMEM_AUTO) then
     memorySampler = ut_sysMemAutoDetect()
  else
     memorySampler = memorySamplerArg
  end if

  do m = 1, AllFns
     c_num_stats_fn(m) = C_NULL_FUNPTR
     c_memory_usage_fn(m) = C_NULL_FUNPTR
  end do

  numFns = 0
  if (iand(memorySampler, UT_SYSMEM_PROC) == UT_SYSMEM_PROC) then
     numFns = numFns + 1
     c_num_stats_fn(numFns) = c_funloc(ut_sysMemProcNumStats)
     c_memory_usage_fn(numFns) = c_funloc(ut_sysMemProc)
  end if
  if (iand(memorySampler, UT_SYSMEM_RUSAGE) == UT_SYSMEM_RUSAGE) then
     numFns = numFns + 1
     c_num_stats_fn(numFns) = c_funloc(ut_sysMemRusageNumStats)
     c_memory_usage_fn(numFns) = c_funloc(ut_sysMemRusage)
  end if
  if (iand(memorySampler, UT_SYSMEM_MALLINFO) == UT_SYSMEM_MALLINFO) then
     numFns = numFns + 1
     c_num_stats_fn(numFns) = c_funloc(ut_sysMemMallinfoNumStats)
     c_memory_usage_fn(numFns) = c_funloc(ut_sysMemMallinfo)
  end if
  if (iand(memorySampler, UT_SYSMEM_BGKERNEL) == UT_SYSMEM_BGKERNEL) then
     numFns = numFns + 1
     c_num_stats_fn(numFns) = c_funloc(ut_sysMemBGKernelNumStats)
     c_memory_usage_fn(numFns) = c_funloc(ut_sysMemBGKernel)
  end if


  numStats = 0
  indexPtr = 1
  do m = 1, numFns
     if ( c_associated(c_memory_usage_fn(m)) .and. &
          c_associated(c_num_stats_fn(m)) ) then

        !Find the number of statistics that the memory sampler will return.
        call c_f_procpointer(c_num_stats_fn(m), num_stats_fn)
        fnStats = num_stats_fn(verbosity)

        if (numStats + fnStats > meminfoSize) then
           !Not enough space in our memory data structure to hold
           !any more statistics.  Exit the do loop.
           exit
        else
           !There is sufficient space.  Nullify the description field for each
           !element because we will soon test if it is null.
           numStats = numStats + fnStats
           do i = indexPtr, numStats
              c_meminfo(i) % description = C_NULL_PTR
           end do
        end if

        !Store the memory statistics in c_meminfo
        call c_f_procpointer(c_memory_usage_fn(m), memory_usage_fn)
        call memory_usage_fn(c_meminfo(indexPtr), fnStats, verbosity)
        indexPtr = indexPtr + fnStats
     end if
  end do


  do i = 1, numStats
     if ( c_associated(c_meminfo(i) % description) ) then
        call c_f_pointer(c_meminfo(i) % description, msgPtr, (/MAX_STRING_LENGTH/))
        meminfo(i) % description = ""
        do j = 1, MAX_STRING_LENGTH
           c = msgPtr(j)
           if (c == C_NULL_CHAR) then
              exit
           else
              meminfo(i) % description(j:j) = c
           end if
        end do
        meminfo(i) % measurement = c_meminfo(i) % measurement
     else
        meminfo(i) % description = "meminfo error"
        meminfo(i) % measurement = -1.0
     end if
  end do
end subroutine ut_sysMemStats
