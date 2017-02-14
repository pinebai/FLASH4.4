!!****f* source/IO/IO_startProtonWrite
!!
!! NAME
!!    IO_startProtonWrite
!!
!! SYNOPSIS
!!
!!    IO_startProtonWrite()
!!
!! DESCRIPTION
!!
!!   This routine reopens the plot file so that proton data can be
!!   written to it. It also creates the extendible 'ProtonData'
!!   dataset in the HDF5 file by calling io_h5create_dataset_protons.
!!
!!***

subroutine IO_startProtonWrite ()

  implicit none

end subroutine IO_startProtonWrite
