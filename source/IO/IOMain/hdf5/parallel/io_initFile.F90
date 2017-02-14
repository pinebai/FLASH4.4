!!****if* source/IO/IOMain/hdf5/parallel/io_initFile
!!
!! NAME
!!  io_initFile
!!
!! SYNOPSIS
!!
!!  io_initFile(integer(in)      :: filenum,
!!              integer(out)     :: fileID, 
!!              character(out)   :: filename(MAX_STRING_LENGTH),
!!              character(in)    :: outputType
!!              logical(in)      :: forced)
!!
!! DESCRIPTION
!!
!!  
!!  Initialized the hdf5 file
!!
!! ARGUMENTS
!!
!!
!!  filenum - number order of output file, 1,2,768 
!!
!!  fileID - file handle returned from hdf5 init file calls
!!
!!  filename - name of the returned file made up of the outputType, filenum and basename
!!
!!  outputType - string indicating output type, usually "chk" or "plt_cnt" or "part"
!!
!!  forced - .true. if this file is "forced."
!!
!!***


subroutine io_initFile( filenum, fileID, filename, outputType, forced)

  use IO_data, ONLY : io_comm, io_outputSplitNum
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none
  
#include "constants.h"

  integer, intent(IN) :: filenum
  integer, intent(INOUT) :: fileID
  character(len=*), intent(in) :: outputType
  character (len=MAX_STRING_LENGTH), intent(inout) :: filename
  logical, intent(IN) :: forced
  integer :: existing

  call io_getOutputName(filenum, "hdf5", outputType, filename, forced)

  existing = 0
  fileID = -1
  call io_h5init_file(fileID, filename, io_comm, io_outputSplitNum, existing)
  if(fileID == -1) then
     print *, "Error: unable to initialize file"
     call Driver_abortFlash("unable to initialize hdf5 file")
  end if

end subroutine io_initFile
