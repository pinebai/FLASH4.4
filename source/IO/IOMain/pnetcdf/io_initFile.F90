!!****if* source/IO/IOMain/pnetcdf/io_initFile
!!
!! NAME
!!  io_initFile
!!
!! SYNOPSIS
!!
!!  io_initFile(integer(in)      :: filenum,
!!              integer(out)     :: fileID, 
!!              character(out)   :: filename(MAX_STRING_LENGTH),
!!              character(in)    :: outputType,
!!              logical(in)      :: forced)
!!
!! DESCRIPTION
!!
!!  
!!  Initialized the pnetcdf file
!!
!! ARGUMENTS
!!
!!  filenum - number order of output file, 1,2,768 
!!
!!  fileID - file handle returned from hdf5 init file calls
!!
!!  filename - name of the returned file made up of the outputType, filenum and basename
!!
!!  outputType - string indicating output type, usually "chk" or "plt_cnt" or "part"
!!
!!  forced - .true. if file is considered "forced."
!!***


subroutine io_initFile( filenum, fileID, filename, outputType, forced)

  use IO_data, ONLY : io_comm, io_outputSplitNum

  implicit none
  
#include "constants.h"

  character (len=*), intent(in) :: outputType
  integer, intent(IN) :: fileID, filenum
  character (len=MAX_STRING_LENGTH), intent(inout) :: filename
  logical, intent(IN) :: forced
  
  call io_getOutputName(filenum, "ncmpi", outputType, filename, forced)

  call io_ncmpi_initialize_file(fileID, filename, io_comm, io_outputSplitNum)



end subroutine io_initFile
