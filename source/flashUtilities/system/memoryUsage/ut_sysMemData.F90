!!****if* source/flashUtilities/system/memoryUsage/ut_sysMemData
!!
!! NAME
!!  ut_sysMemData
!!
!! SYNOPSIS
!!  use ut_sysMemData
!!   
!!***

#include "constants.h"

module ut_sysMemData
  implicit none
  type meminfo_t
     real :: measurement
     character(len=MAX_STRING_LENGTH) :: description
  end type meminfo_t

  type memsummary_t
     real :: measurement, min, max, avg
     integer :: min_rank, max_rank
     character(len=MAX_STRING_LENGTH) :: description
  end type memsummary_t
end module ut_sysMemData
