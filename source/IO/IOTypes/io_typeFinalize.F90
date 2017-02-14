!!****if* source/IO/IOTypes/io_typeFinalize
!!
!! NAME
!!  io_typeFinalize
!!
!! SYNOPSIS
!!
!!  io_typeFinalize()
!!
!! DESCRIPTION
!!
!!  Cleans up IOTypes sub unit.
!!
!! ARGUMENTS
!!
!!***

#include "constants.h"
#include "Flash.h"

subroutine io_typeFinalize()
#ifdef USE_IO_C_INTERFACE
  use io_c_type_interface, ONLY : io_free_grid_mpi_types
#ifdef FLASH_IO_PNETCDF
  use io_c_type_interface, ONLY : io_ncmpi_nonblocking_finalize
#endif
#endif
  implicit none

  call io_free_grid_mpi_types()

#ifdef FLASH_IO_PNETCDF
  call io_ncmpi_nonblocking_finalize()
#endif
end subroutine io_typeFinalize
