!!****if* source/IO/IOMain/hdf5/IO_endProtonWrite
!!
!! NAME
!!
!!  IO_endProtonWrite
!!
!! SYNOPSIS
!!
!!  call IO_endProtonWrite ()
!!
!! DESCRIPTION
!!
!!   This subroutine is called after all of the proton data has been written
!!   to the HDF5 plot file. It simply closes the plot file.
!!
!!***

subroutine IO_endProtonWrite ()

  use IO_data,          ONLY: io_protonFileID
  use Driver_interface, ONLY: Driver_abortFlash

  implicit none
!
!
!    ...Close the HDF5 plot file.
!
!
  call io_h5close_file (io_protonFileID)

end subroutine IO_endProtonWrite
