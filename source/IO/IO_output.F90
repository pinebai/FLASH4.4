!!****f* source/IO/IO_output
!!
!! NAME
!!
!!  IO_output
!!
!!
!! SYNOPSIS
!!
!!
!!  IO_output(real(in)              :: simTime, 
!!            real(in)              :: dt,
!!            integer(in)           :: nstep,
!!            integer(in)           :: nbegin,
!!            logical(out)          :: endRun,
!!            integer(in), optional :: outputType)
!!
!!
!! DESCRIPTION
!!
!!  This routine is called after every sweep to check if it is time to
!!  output.  Checkpoints and plotfiles and .dat files
!!  are handled through this function.
!!  
!!  
!!  A checkpoint can be triggered in a few different ways.  First, enough
!!  wall clock time has elapsed since the code has started running that
!!  we want to force a checkpoint now.  This is controlled by the
!!  wall_clock_checkpoint runtime parameter.  This is always timed from
!!  the start of execution of the code, and not since the last checkpoint,
!!  so we can ensure that we get a checkpoint dumped right before the
!!  queue window closes.  Checkpoints are also produced after a given
!!  amount of simulation time has elapsed -- this is controlled by the
!!  checkpointFileIntervalTime runtime parameter.  Finally, a checkpoint can be forced in
!!  to be produced every checkpointFileIntervalStep timesteps or by creating .dump_checkpoint
!!  file in the execution directory.
!!
!!  Plotfiles are produced equally spaced in simulation time, tplot time
!!  units apart. They can also be produced on demand by temporarily creating
!!  .dump_plotfile execution control file.
!!  
!!
!!  After every sweep, IO_writeIntegralQuantities is called to compute some global
!!  quantities and write them to the flash.dat file.
!!
!!  A checkpoint is given a name according to the basename (specified
!!  via the basenm runtime parameter) and the filetype (pnetcdf,
!!  HDF5).  Each separate checkpoint is given a unique number suffix.  The
!!  value to start with (or restart from if this simulation was a restart
!!  from a previous checkpoint) is set by the checkpointFileNumber runtime parameter.
!!  Or if it is a restart, checkpointFileNumber can be saved as a scalar in the checkpoint
!!  file.
!!
!!  Execution will be aborted immediately if .kill file is present
!!  in the execution directory. If .dump_restart is present then a
!!  checkpoint file is saved before aborting execution.
!!
!!  For more detail about runtime parameters controlling IO output
!!  look at the IO/common Config file or the setup_params text file
!!  written to the object directory.
!!
!! ARGUMENTS
!!
!!  simTime - simulation time
!!  dt - timestep
!!  nstep - current time step number
!!  nbegin - beginning time step number, if starting from scratch this is 0,
!!      it could be different in case of restart
!!  endRun - will be set to .TRUE. on return if existence of a .dump_restart or
!!      a .kill file was detected; .FALSE. otherwise.
!!  outputType - an integer that denotes the type of output files you expect to 
!!               generate with this call. If this argument is omitted IO_output
!!               will test to see if it is time to output every kind of file.
!!
!! NOTES
!!
!!  Variables that start with "io_", like io_checkpointFileNumber and io_checkpointFileIntervalTime
!!  are located in the IO_data fortran 
!!  module.  The "io_" is meant to indicated that this variable has
!!  IO Unit scope.  Variable without the "io_" in front are local
!!  variables to this subroutine.
!!
!!  The outputType argument has a series of constants declared in constants.h which
!!  allow you to state which filed you wish to have output (povided that the normal 
!!  output conditions are met) by this call.  The options are:
!!
!!  CHECKPOINT_FILE_ONLY 
!!  PLOTFILE_ONLY 
!!  PRTICLE_FILE_ONLY 
!!  CHECKPOINT_AND_PLOTFILE 
!!  CHECKPOINT_AND_PARTICLEFILE 
!!  PLOTFILE_AND_PARTICLEFILE 
!!  ALL_FILES 
!!
!!  If the optional outputType arguemnt is omitted, ALL_FILES is assumed.
!!
!! SIDE EFFECTS
!!
!!  For visualization purposes, the data is restricted up the entire tree,
!!  so the data is valid on all levels.  This allows multi-resolution vis
!!  techniques to be applied.  Immediately after checkpointing or writing
!!  a plotfile is the only time the data is guaranteed to be valid on all
!!  levels.
!!
!! SEE ALSO
!!   IO_writeCheckpoint, IO_writeIntegralQuantities, IO_writeParticles,
!!   IO_writePlotfile, IO_writeUserArray, IO_setScalar
!!
!!***

!!  Variables that start with "io_", like io_checkpointFileNumber and
!!  io_checkpointFileIntervalTime are located in the IO_data fortran
!!  module.  The "io_" is meant to indicated that this variable has IO
!!  Unit scope.  Variable without the "io_" in front are local
!!  variables to this subroutine


subroutine IO_output(simTime, dt, nstep, nbegin, endRun, outputType)

  use IO_interface, ONLY : IO_writeIntegralQuantities

  implicit none

  real, intent(in) :: simTime, dt
  integer, intent(in) :: nstep, nbegin
  logical, intent(out) :: endRun
  integer, intent(in), optional :: outputType

  call IO_writeIntegralQuantities( 0, simTime)

  endRun = .FALSE.

end subroutine IO_output
