!!****ih* source/flashUtilities/system/memoryUsage/f2003/ut_sysMemCInterface
!!
!! NAME
!!  ut_sysMemCInterface 
!!
!! SYNOPSIS
!!  use ut_sysMemCInterface
!!
!!***

#include "constants.h"

module ut_sysMemCInterface
  implicit none

  interface
     function ut_sysMemAutoDetect() &
          bind(c,name='ut_sysMemAutoDetect')
       use iso_c_binding, ONLY : c_int
       integer(c_int) :: ut_sysMemAutoDetect
     end function ut_sysMemAutoDetect
  end interface

  abstract interface
     subroutine ut_sysMem &
          (meminfo,meminfoSize,verbosity) &
          bind(c)
       use iso_c_binding, ONLY : c_int
       use ut_sysMemCData, ONLY : c_meminfo_t
       type(c_meminfo_t), dimension(*), intent(inout) :: meminfo
       integer(c_int), value, intent(in) :: meminfoSize, verbosity
     end subroutine ut_sysMem
  end interface

  abstract interface
     function ut_sysMemNumStats (verbosity) &
          bind(c)
       use iso_c_binding, ONLY : c_int
       integer(c_int), value, intent(in) :: verbosity
       integer(c_int) :: ut_sysMemNumStats
     end function ut_sysMemNumStats
  end interface

  interface
     subroutine ut_sysMemMallinfo &
          (meminfo,meminfoSize,verbosity) &
       bind(c,name='ut_sysMemMallinfo')
       use iso_c_binding, ONLY : c_int
       use ut_sysMemCData, ONLY : c_meminfo_t
       type(c_meminfo_t), dimension(*), intent(inout) :: meminfo
       integer(c_int), value, intent(in) :: meminfoSize, verbosity
     end subroutine ut_sysMemMallinfo
  end interface

  interface
     function ut_sysMemMallinfoNumStats (verbosity) &
          bind(c,name='ut_sysMemMallinfoNumStats')
       use iso_c_binding, ONLY : c_int
       integer(c_int), value, intent(in) :: verbosity
       integer(c_int) :: ut_sysMemMallinfoNumStats
     end function ut_sysMemMallinfoNumStats
  end interface

  interface
     subroutine ut_sysMemProc &
          (meminfo,meminfoSize,verbosity) &
          bind(c,name='ut_sysMemProc')
       use iso_c_binding, ONLY : c_int
       use ut_sysMemCData, ONLY : c_meminfo_t
       type(c_meminfo_t), dimension(*), intent(inout) :: meminfo
       integer(c_int), value, intent(in) :: meminfoSize, verbosity
     end subroutine ut_sysMemProc
  end interface

  interface
     function ut_sysMemProcNumStats (verbosity) &
          bind(c,name='ut_sysMemProcNumStats')
       use iso_c_binding, ONLY : c_int
       integer(c_int), value, intent(in) :: verbosity
       integer(c_int) :: ut_sysMemProcNumStats
     end function ut_sysMemProcNumStats
  end interface

  interface
     subroutine ut_sysMemBGKernel &
          (meminfo,meminfoSize,verbosity) &
          bind(c,name='ut_sysMemBGKernel')
       use iso_c_binding, ONLY : c_int
       use ut_sysMemCData, ONLY : c_meminfo_t
       type(c_meminfo_t), dimension(*), intent(inout) :: meminfo
       integer(c_int), value, intent(in) :: meminfoSize, verbosity
     end subroutine ut_sysMemBGKernel
  end interface

  interface
     function ut_sysMemBGKernelNumStats (verbosity) &
          bind(c,name='ut_sysMemBGKernelNumStats')
       use iso_c_binding, ONLY : c_int
       integer(c_int), value, intent(in) :: verbosity
       integer(c_int) :: ut_sysMemBGKernelNumStats
     end function ut_sysMemBGKernelNumStats
  end interface

  interface
     subroutine ut_sysMemRusage &
          (meminfo,meminfoSize,verbosity) &
          bind(c,name='ut_sysMemRusage')
       use iso_c_binding, ONLY : c_int
       use ut_sysMemCData, ONLY : c_meminfo_t
       type(c_meminfo_t), dimension(*), intent(inout) :: meminfo
       integer(c_int), value, intent(in) :: meminfoSize, verbosity
     end subroutine ut_sysMemRusage
  end interface

  interface
     function ut_sysMemRusageNumStats (verbosity) &
          bind(c,name='ut_sysMemRusageNumStats')
       use iso_c_binding, ONLY : c_int
       integer(c_int), value, intent(in) :: verbosity
       integer(c_int) :: ut_sysMemRusageNumStats
     end function ut_sysMemRusageNumStats
  end interface

end module ut_sysMemCInterface
