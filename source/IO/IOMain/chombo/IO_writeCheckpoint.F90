!!****if* source/IO/IOMain/chombo/IO_writeCheckpoint
!!
!! NAME
!!
!!  IO_writeCheckpoint
!!
!!
!! SYNOPSIS
!!
!!  IO_writeCheckpoint()
!!
!!
!!
!! DESCRIPTION
!!
!!  This is a generic call to write the important simulation data to a
!!  checkpoint file.  A checkpoint file writes a few different types of
!!  data to a file, first the physical data like, temperature, pressure, density
!!  etc. which are stored in each cell on the grid.  
!!  Second, in order to recreate the simulation from a checkpoint file a
!!  number of other single quantities are needed as well.  We call these
!!  scalar values which include simTime, dt, nstep, globalNumBlocks etc.
!!  We also store descriptive strings that describe the simulation run.
!!
!!  The same IO_writeCheckpoint routine is called regardless of the type of
!!  file being written, (such as hdf5 parallel, hdf5 serial or pnetcdf)
!!  IO_writeCheckpoint prepares the Grid_ioData (like getting the
!!  globalNumBlocks) and collects the scalars wanting to be checkpointed
!!  from each unit. 
!!  IO_writeCheckpoint then calls four methods, io_initFile, io_writeData, 
!!  IO_writeParticles and io_closeFile.  Each of these routines _is_ specific
!!  to the type of io library used and have their own implementation.  
!!  In addition, io_writeData has its own
!!  implementation for io library and type of grid (UG, Paramesh, or other)
!!
!!
!!  In FLASH IO_writeCheckpoint is called from IO_output (or IO_outputInitial or
!!  IO_outputFinal) IO_output checks to see if enough wall clock time,
!!  simTim, or nsteps has passed to checkpoint.
!!
!!  We have put IO_writeCheckpoint in the API because a user may want to write
!!  a checkpoint at another time or for another reason without having to go through
!!  IO_output.  For most flash users IO_writeCheckpoint will only ever be
!!  called through IO_output
!!
!! ARGUMENTS
!! 
!!
!! NOTES
!!
!!  For those familiar with FLASH2, breaking up the checkpoint routine into
!!  these four different methods is a change.  Because FLASH3 now supports
!!  different grid packages and we are committed to supporting both
!!  hdf5 and parallel netCDF having each grid and io library writing its
!!  own checkpoint file proved to be a lot of code duplication.  We believe
!!  that while dividing up the checkpoint routines created more files it 
!!  will in the end be easier to maintain.
!!
!!***
subroutine IO_writeCheckpoint()
  use Grid_interface, ONLY : Grid_dump
  use IO_interface, ONLY : IO_updateScalars
  use Logfile_interface, ONLY : Logfile_stamp
  use iso_c_binding!, ONLY : C_NULL_CHAR, c_loc, c_double, c_int, c_char, c_bool
  use flash_ftypes
  use chombo_f_c_interface, ONLY : ch_write_checkpoint
  use Driver_interface, only: Driver_getSimTime, Driver_getDt
  use IO_data
  implicit none

#include "constants.h"  
#include "Flash.h"  

  character (len=MAX_STRING_LENGTH) :: filename
  ! for logfile output
  character(len=MAX_STRING_LENGTH), dimension(1,2) :: strBuff
  real(c_double) :: simTime, dt
  type(named_vals_t) :: sc_vals, rp_vals
  integer :: i
  
  call io_getOutputName(io_checkpointFileNumber, "chombo_hdf5", "_chk_", filename, .false.)

  if (io_globalMe == MASTER_PE) then
     write (strBuff(1,1), "(A)") "type"
     write (strBuff(1,2), "(A)") "checkpoint"
     call Logfile_stamp(strBuff, 1, 2, "[IO_writeCheckpoint] close")
     print *, '*** Wrote checkpoint file ***' 
  end if

  call Driver_getSimTime(simTime)
  call Driver_getDt(dt)
  
  call IO_updateScalars()
  call io_prepareSimInfo()
  call io_prepareListsWrite()
  
  ! runtime parameters
  rp_vals%real_count = io_numRealParms
  rp_vals%real_names = c_loc(io_realParmNames)
  rp_vals%real_vals = c_loc(io_realParmValues)
  
  rp_vals%int_count = io_numIntParms
  rp_vals%int_names = c_loc(io_intParmNames)
  rp_vals%int_vals = c_loc(io_intParmValues)
  
  rp_vals%str_count = io_numStrParms
  rp_vals%str_names = c_loc(io_strParmNames)
  rp_vals%str_vals = c_loc(io_strParmValues)
  
  rp_vals%log_count = io_numLogParms
  rp_vals%log_names = c_loc(io_logParmNames)
  rp_vals%log_vals = c_loc(io_logToIntParmValues)
  
  ! scalars
  sc_vals%real_count = io_numRealScalars
  sc_vals%real_names = c_loc(io_realScalarNames)
  sc_vals%real_vals = c_loc(io_realScalarValues)
  
  sc_vals%int_count = io_numIntScalars
  sc_vals%int_names = c_loc(io_intScalarNames)
  sc_vals%int_vals = c_loc(io_intScalarValues)
  
  sc_vals%str_count = io_numStrScalars
  sc_vals%str_names = c_loc(io_strScalarNames)
  sc_vals%str_vals = c_loc(io_strScalarValues)
  
  sc_vals%log_count = io_numLogScalars
  sc_vals%log_names = c_loc(io_logScalarNames)
  sc_vals%log_vals = c_loc(io_logToIntScalarValues)
  
  call ch_write_checkpoint(trim(filename)//C_NULL_CHAR, simTime, dt, sc_vals, rp_vals)
  
  call io_finalizeListsWrite()
  io_checkpointFileNumber = io_checkpointFileNumber + 1

end subroutine IO_writeCheckpoint
