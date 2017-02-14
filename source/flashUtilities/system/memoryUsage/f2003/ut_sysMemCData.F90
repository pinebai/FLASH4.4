!!****if* source/flashUtilities/system/memoryUsage/f2003/ut_sysMemCData
!!
!! NAME
!!  ut_sysMemCData
!!
!! SYNOPSIS
!!  use ut_sysMemCData
!!   
!!***

#include "constants.h"

module ut_sysMemCData
  use iso_c_binding, ONLY : c_double, c_ptr
  implicit none

  type, bind(c) :: c_meminfo_t
     real(c_double) :: measurement
     type(c_ptr) :: description
  end type c_meminfo_t
end module ut_sysMemCData
