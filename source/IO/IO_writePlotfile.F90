!!****f* source/IO/IO_writePlotfile
!!
!!
!! NAME
!!
!!  IO_writePlotfile
!!
!!
!! SYNOPSIS
!!
!!  IO_writePlotfile(logical(in), optional :: forced)
!!
!!
!!
!! DESCRIPTION
!!
!!  This is a generic call to write the important simulation data to a
!!  plotfile file.  A plotfile file writes a few different types of
!!  data to a file, first the physical data like, temperature,
!!  pressure, density etc.  for all cells on the grid.  Secondly, in
!!  order to recreate the simulation from a plotfile file a number of
!!  other single quantities are needed as well.  We call these scalar
!!  values which include simTime, dt, nstep, globalNumBlocks etc.  We
!!  also store descriptive strings that describe the simulation run.
!!
!!  The same IO_writePlotfile routine is called regardless of the type
!!  of file being written, (such as hdf5 parallel, hdf5 serial or
!!  pnetcdf) IO_writePlotfile prepares the Grid_ioData (like getting
!!  the globalNumBlocks) and collects the scalars wanting to be stored
!!  from each unit. IO_writePlotfile then calls four methods,
!!  io_initFile, io_writeData, and io_closeFile.  Each of these
!!  routines _is_ specific to the type of io library used and have
!!  their own implementation.  In addition, io_writeData has its own
!!  implementation for io library and type of grid (UG, Paramesh, or
!!  other)
!!
!!  Since plotfiles are used for visualization purposes, and not for
!!  restarting a run, to keep plotfile sizes manageable, data is
!!  written out in single precision
!!
!!  In FLASH IO_writePlotfile is called from IO_output (or
!!  IO_outputInitial or IO_outputFinal) IO_output checks whether it is
!!  time to output a plotfile.  The runtime parameters that control
!!  writing plotfiles are tplot, the simulation time between plotfiles
!!  and nplot, the number of timesteps between plotfiles.
!!  
!!
!!  We have put IO_writePlotfile in the API because a user may want to
!!  write a plotfile at another time or for another reason without
!!  having to go through IO_output.  For most flash users
!!  IO_writePlotfile will only ever be called through IO_output.
!!
!! ARGUMENTS
!! 
!!  forced - .true. if this is a "forced" plotfile.  Default is .false.
!!
!! NOTES
!!  
!!  We have added functionality to compress plotfile data further.
!!  Byte packing is currently only available with hdf5 and the uniform
!!  grid.  To use byte packing set the runtime parameter bytePack to
!!  .true. in your flash.par
!!
!!  For those familiar with FLASH2, breaking up the plotfile routine into
!!  these four different methods is a change.  Because FLASH3 now supports
!!  different grid packages and we are committed to supporting both
!!  hdf5 and parallel netCDF having each grid and io library writing its
!!  own plotfile file proved to be a lot of code duplication.  We believe
!!  that while dividing up the plotfile routines created more files it 
!!  will in the end be easier to maintain.
!!
!! SEE ALSO
!!  IO_output
!! 
!!
!!
!!***


subroutine IO_writePlotfile( forced)

  implicit none

  logical, intent(in), optional :: forced

  return
end subroutine IO_writePlotfile




