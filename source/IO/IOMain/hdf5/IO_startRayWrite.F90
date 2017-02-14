!!****if* source/IO/IOMain/hdf5/IO_startRayWrite
!!
!! NAME
!!    IO_startRayWrite
!!
!! SYNOPSIS
!!
!!    IO_startRayWrite()
!!
!! DESCRIPTION
!!
!!   This routine reopens the plot file so that laser rays can be
!!   written to it. It also creates the extendible RayData dataset in
!!   the HDF5 file by calling io_h5create_raydset.
!!
!!***
subroutine IO_startRayWrite()
  use IO_data, ONLY: io_wrotePlot,       &
                     io_oldPlotFileName, &
                     io_meshComm,        &
                     io_outputSplitNum,  &
                     io_rayFileID

  use Driver_interface, ONLY: Driver_abortFlash
  implicit none

#include "constants.h"

  integer :: existing

  if(.not. io_wrotePlot) then
     call Driver_abortFlash("[IO_startRayWrite] Rays can only be written after a plot")
  end if

  ! Re-open the HDF5 plot file:
  existing = 1
  io_rayFileID = -1
  call io_h5init_file(io_rayFileID, io_oldPlotFileName, io_meshComm, io_outputSplitNum, existing)
  if(io_rayFileID == -1) then
     call Driver_abortFlash("[IO_writeRays] unable to open hdf5 file: " // &
          trim(io_oldPlotFileName))
  end if

  ! Create an extendible dataset to store ray data:
  call io_h5create_raydset(io_rayFileID)          

end subroutine IO_startRayWrite
