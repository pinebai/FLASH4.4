!!****if* source/IO/IOMain/hdf5/serial/io_closeFile
!!
!! NAME
!!  io_closeFile
!!
!! SYNOPSIS
!!
!!  io_closeFile(integer(in) :: fileID)
!!
!! DESCRIPTION
!!  
!!  closes the hdf5 file      
!!
!! ARGUMENTS
!!
!!  fileID - integer file identifier for hdf5 file
!!
!!***

subroutine io_closeFile( fileID)

  use IO_data, ONLY : io_globalMe
  implicit none

#include "constants.h"

  integer, intent(in) :: fileID

  if(io_globalMe /= MASTER_PE) then
      return
  end if

  call io_h5close_file(fileID)

end subroutine io_closeFile
