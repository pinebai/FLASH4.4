!!****if* source/IO/IOMain/hdf5/IO_startProtonWrite
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

  use IO_data,          ONLY: io_wrotePlot,       &
                              io_oldPlotFileName, &
                              io_meshComm,        &
                              io_outputSplitNum,  &
                              io_protonFileID

  use Driver_interface, ONLY: Driver_abortFlash

  implicit none

#include "constants.h"

  integer :: existing

  if (.not. io_wrotePlot) then
       call Driver_abortFlash("[IO_startProtonWrite] Protons can only be written after a plot")
  end if
!
!
!    ...Re-open the HDF5 plot file.
!
!
  existing = 1              ! tells the init file routine that the plot file already exists
  io_protonFileID = -1      ! overwritten by the file ID returned from the init file routine

  call io_h5init_file (io_protonFileID,    &
                       io_oldPlotFileName, &
                       io_meshComm,        &
                       io_outputSplitNum,  &
                       existing            )  ! -1 means 'open' the file, not 'create'

  if (io_protonFileID == -1) then
      call Driver_abortFlash("[IO_startProtonWrite] unable to open hdf5 file: " // &
                              trim (io_oldPlotFileName))
  end if
!
!
!    ...Create an extendible 'ProtonData' dataset inside the HDF5 plot file,
!       to be able to store proton data.
!
!
  call io_h5create_dataset_protons (io_protonFileID)          

end subroutine IO_startProtonWrite
