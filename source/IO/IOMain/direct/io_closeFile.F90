!!****if* source/IO/IOMain/direct/io_closeFile
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
!!  fileID - integer file identifier
!!
!!***

subroutine io_closeFile(fileID)

  implicit none

  integer, intent(in) :: fileID

  close(fileID)

end subroutine io_closeFile
