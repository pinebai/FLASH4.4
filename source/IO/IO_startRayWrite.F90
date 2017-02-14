!!****f* source/IO/IO_startRayWrite
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
  implicit none

end subroutine IO_startRayWrite
