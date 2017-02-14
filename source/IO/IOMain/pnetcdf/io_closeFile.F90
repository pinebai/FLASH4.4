!!****if* source/IO/IOMain/pnetcdf/io_closeFile
!!
!! NAME
!!
!!  io_closeFile
!!
!! SYNOPSIS
!!
!!  io_closeFile(integer, intent(in)  :: fileid)
!!
!! DESCRIPTION
!!   closes the pnetcdf file
!!
!! ARGUMENTS
!!
!!   fileid : pnetcdf file identifier
!!
!! 
!!
!!
!!***

subroutine io_closeFile(fileID)

  implicit none

  integer, intent(in) :: fileID

  call io_ncmpi_close_file(fileID)

end subroutine io_closeFile
