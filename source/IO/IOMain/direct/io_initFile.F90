!!****if* source/IO/IOMain/direct/io_initFile
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
!!              logical(in)           :: forced
!!
!! DESCRIPTION
!!
!!  
!!  Initialized the basic output file
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
!!  forced - .true. if this is considered a "forced" output file.
!!***


subroutine io_initFile(filenum, fileID, filename, outputType, forced)

  use IO_data, ONLY : io_baseName, io_outputDir, io_restart, io_meshMe

  implicit none
  
#include "constants.h"
#include "Flash.h"

  integer, intent(in) :: filenum
  integer, intent(inout) :: fileID
  character(len=*), intent(in) :: outputType
  character (len=MAX_STRING_LENGTH), intent(inout) :: filename
  character (len=MAX_STRING_LENGTH) :: firstOutputDir, filenameFirst
  logical, intent(in) :: forced

  integer :: ioLun = 28

  ! create a character variable to hold the string representation of the block
  ! number.  Note this is set to be 4 characters long (i.e. max = 9999).  
  character (len=6) ::  fnumStr 
  character (len=6) ::  splitFileNumStr
  integer :: pos
  integer :: myPE
  myPE = io_meshMe

  ! if io_outputDir isn't empty, and it doesn't have a 
  ! directory seperator at the end, add the directory seperator.
  pos = index(io_baseName, ' ')
  write (fnumStr, '(i6.6)') filenum

  pos = index(io_outputDir, ' ')
  if (pos.gt.1) then
     if (io_outputDir(pos-1:pos-1).ne."/") then
        io_outputDir(pos:pos) = "/"
     end if
  end if
  

  ! if io_outputDir is current directory specified 
  ! with a './', just get rid of it.
  if (io_outputDir == "./") then
     io_outputDir = ""
  end if

  filename = " "
  write (fnumStr, '(i6.6)') filenum

  pos = index(io_baseName, ' ')

  !only if we are splitting files do we adjust the filename
  write (splitFileNumStr, '(i6.6)') myPE

  if(forced) then
     filename = io_baseName(:pos-1) //"forced_"// 'direct' // outputType // 's'// splitFileNumStr // '_' // fnumStr 
  else
     filename = io_baseName(:pos-1) // 'direct' // outputType // 's'// splitFileNumStr // '_' // fnumStr 
  end if
  pos = index(io_outputDir, ' ')
  filename = io_outputDir(:pos-1) // filename 


  open(unit=ioLun,file=filename,form="unformatted")

  fileID = ioLun


end subroutine io_initFile
