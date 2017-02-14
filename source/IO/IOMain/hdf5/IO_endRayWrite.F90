!!****if* source/IO/IOMain/hdf5/IO_endRayWrite
!!
!! NAME
!!
!!  IO_endRayWrite
!!
!! SYNOPSIS
!!
!!  call IO_endRayWrite ()
!!
!! DESCRIPTION
!!
!!   This subroutine is called after all of the ray paths have moved
!!   to the plot file. It simply closes the plot file.
!!
!!***
subroutine IO_endRayWrite()
  use IO_data, ONLY: io_wrotePlot,       &
                     io_oldPlotFileName, &
                     io_outputSplitNum,  &
                     io_rayFileID

  use Driver_interface, ONLY: Driver_abortFlash
  implicit none

  if(.not. io_wrotePlot) then
     call Driver_abortFlash("[IO_endRayWrite] Rays can only be written after a plot")
  end if

  ! Close the HDF5 plot file:
  call io_h5close_file(io_rayFileID)

end subroutine IO_endRayWrite
