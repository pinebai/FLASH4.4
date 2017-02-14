!!****if* source/IO/IOMain/chombo/IO_readCheckpoint
!!
!! NAME
!!
!!  IO_readCheckpoint
!!
!!
!! SYNOPSIS
!!
!!  IO_readCheckpoint()
!!
!!
!!
!! DESCRIPTION
!!
!!  IO_readCheckpoint is a generic subroutine that retrieves
!!  the unklabels and then calls io_readData
!!  which is specific to pnetcdf, hdf5 and
!!  the necessary grid, UG, paramesh etc.
!!  io_readData reads a checkpoint file and reinitializes
!!  grid and scalar values to resume the run
!!  
!!
!!***

subroutine IO_readCheckpoint()
  use iso_c_binding
  use chombo_f_c_interface
  use flash_ftypes
  use Driver_interface, ONLY : Driver_abortFlash
  use Driver_data, ONLY : dr_dt, dr_simTime
  use IO_data!, only: io_checkpointFileNumber
  implicit none

#include "constants.h"  
#include "Flash.h"  

  character(kind=c_char, len=MAX_STRING_LENGTH) :: filename
  real(c_double) :: simtime, dt
  type(named_vals_t) :: nv_sc, nv_rp
  real(c_double), pointer :: preal(:)
  integer(c_int), pointer :: pint(:)
  character(kind=c_char, len=MAX_STRING_LENGTH), pointer :: pstr(:)
  
  call io_getOutputName(io_checkpointFileNumber, "chombo_hdf5", "_chk_", filename, .false.)
  
  call ch_read_checkpoint(trim(filename)//C_NULL_CHAR, simtime, dt, nv_sc, nv_rp)
  dr_simTime = simtime
  dr_dt = dt
  
  io_numRealScalars = nv_sc%real_count
  io_numIntScalars = nv_sc%int_count
  io_numStrScalars = nv_sc%str_count
  io_numLogScalars = nv_sc%log_count
  io_numRealParmsPrev = nv_rp%real_count
  io_numIntParmsPrev = nv_rp%int_count
  io_numStrParmsPrev = nv_rp%str_count
  io_numLogParmsPrev = nv_rp%log_count
  call io_prepareListsRead()
  
  ! scalars
  call c_f_pointer(nv_sc%real_names, pstr, (/io_numRealScalars/))
  io_realScalarNames(:) = pstr(:)
  call c_f_pointer(nv_sc%real_vals, preal, (/io_numRealScalars/))
  io_realScalarValues(:) = preal(:)
  call c_free(nv_sc%real_names)
  call c_free(nv_sc%real_vals)
  
  call c_f_pointer(nv_sc%int_names, pstr, (/io_numIntScalars/))
  io_intScalarNames(:) = pstr(:)
  call c_f_pointer(nv_sc%int_vals, pint, (/io_numIntScalars/))
  io_intScalarValues(:) = pint(:)
  call c_free(nv_sc%int_names)
  call c_free(nv_sc%int_vals)
  
  call c_f_pointer(nv_sc%str_names, pstr, (/io_numStrScalars/))
  io_strScalarNames(:) = pstr(:)
  call c_f_pointer(nv_sc%str_vals, pstr, (/io_numStrScalars/))
  io_strScalarValues(:) = pstr(:)
  call c_free(nv_sc%str_names)
  call c_free(nv_sc%str_vals)
  
  call c_f_pointer(nv_sc%log_names, pstr, (/io_numLogScalars/))
  io_logScalarNames(:) = pstr(:)
  call c_f_pointer(nv_sc%log_vals, pint, (/io_numLogScalars/))
  io_logToIntScalarValues(:) = pint(:)
  call c_free(nv_sc%log_names)
  call c_free(nv_sc%log_vals)

  ! runparms
  call c_f_pointer(nv_rp%real_names, pstr, (/io_numRealParmsPrev/))
  io_realParmNamesPrev(:) = pstr(:)
  call c_f_pointer(nv_rp%real_vals, preal, (/io_numRealParmsPrev/))
  io_realParmValuesPrev(:) = preal(:)
  call c_free(nv_rp%real_names)
  call c_free(nv_rp%real_vals)
  
  call c_f_pointer(nv_rp%int_names, pstr, (/io_numIntParmsPrev/))
  io_intParmNamesPrev(:) = pstr(:)
  call c_f_pointer(nv_rp%int_vals, pint, (/io_numIntParmsPrev/))
  io_intParmValuesPrev(:) = pint(:)
  call c_free(nv_rp%int_names)
  call c_free(nv_rp%int_vals)
  
  call c_f_pointer(nv_rp%str_names, pstr, (/io_numStrParmsPrev/))
  io_strParmNamesPrev(:) = pstr(:)
  call c_f_pointer(nv_rp%str_vals, pstr, (/io_numStrParmsPrev/))
  io_strParmValuesPrev(:) = pstr(:)
  call c_free(nv_rp%str_names)
  call c_free(nv_rp%str_vals)
  
  call c_f_pointer(nv_rp%log_names, pstr, (/io_numLogParmsPrev/))
  io_logParmNamesPrev(:) = pstr(:)
  call c_f_pointer(nv_rp%log_vals, pint, (/io_numLogParmsPrev/))
  io_logToIntParmValuesPrev(:) = pint(:)
  call c_free(nv_rp%log_names)
  call c_free(nv_rp%log_vals)
  
  call io_finalizeListsRead()
  
  io_checkpointFileNumber = io_checkpointFileNumber + 1
  
end subroutine IO_readCheckpoint

