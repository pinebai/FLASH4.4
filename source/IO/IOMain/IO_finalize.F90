!!****if* source/IO/IOMain/IO_finalize
!!
!! NAME
!!  IO_finalize
!!
!! SYNOPSIS
!!
!!  IO_finalize()
!!
!! DESCRIPTION
!!
!!  This function cleans up the IO unit, deallocates memory, etc.
!!
!! ARGUMENTS
!!
!!***

#include "constants.h"
#include "Flash.h"

subroutine IO_finalize()
  implicit none

#ifdef FLASH_IO_EXPERIMENTAL
  call io_typeFinalize()
#endif
end subroutine IO_finalize
