!!****if* source/IO/IOMain/hdf5/io_createDatasets
!!
!! NAME
!!
!!  io_createDatasets
!!
!!
!! SYNOPSIS
!!
!!  io_createDatasets(integer(in) :: fileID,
!!                    integer(in) :: numFileBlks,
!!                    integer(in) :: presentDims) 
!!          
!!
!!
!!
!! DESCRIPTION
!!
!!  Creates the datasets in the HDF5 file.
!!
!!
!! ARGUMENTS
!! 
!!  fileID - integer file identifier for hdf5 file
!!  numFileBlks - number of blocks in the hdf5 file
!!  presentDims - number of dimensions in certain FLASH datasets
!!                (provided for backwards compatibility)
!!
!!***

#include "constants.h"
#include "Flash.h"
#include "io_flash.h"

subroutine io_createDatasets(fileID, numFileBlks, presentDims) 

#ifdef USE_IO_C_INTERFACE
  use iso_c_binding, ONLY : c_char
  use io_c_interface, ONLY : io_create_dataset
#endif
  use IO_data, ONLY : io_globalMe, io_doublePrecision, io_plotfileMetadataDP
#ifdef FLASH_GRID_PARAMESH
  use tree, ONLY : nfaces, nchild
#ifdef FLASH_GRID_PARAMESH3OR4
  use tree, ONLY : MFLAGS
#endif
#endif
  use Timers_interface, ONLY : Timers_start, Timers_stop

  implicit none
  integer, intent(in) :: fileID, numFileBlks, presentDims

  integer :: treeFP
  integer, parameter :: libType = IO_FILE_HDF5
#ifndef USE_IO_C_INTERFACE
  integer, parameter :: c_char = KIND('A')
#endif
#ifndef FLASH_GRID_PARAMESH
  integer, parameter :: nfaces = 2*NDIM, nchild = 2**NDIM
#endif
  character(kind=c_char,len=6) :: bflags_str = "bflags"
  character(kind=c_char,len=3) :: gr_gid_str = "gid"
  character(kind=c_char,len=11) :: coord_str = "coordinates"
  character(kind=c_char,len=9) :: nodetype_str = "node type"
  character(kind=c_char,len=12) :: lrefine_str = "refine level"
  character(kind=c_char,len=11) :: which_child_str = "which child"
  character(kind=c_char,len=16) :: procnumber_str = "processor number"
  character(kind=c_char,len=12) :: bnd_box_str = "bounding box"
  character(kind=c_char,len=10) :: bsize_str = "block size"
  character(kind=c_char,len=10) :: gr_gsurr_blks_str = "gsurr_blks"


  !Floating point type for floating point tree datasets.
  if (io_doublePrecision .or. io_plotfileMetadataDP) then
     treeFP = IO_FLASH_DOUBLE
  else
     treeFP = IO_FLASH_FLOAT
  end if

  call io_create_dataset(io_globalMe, fileID, libType, &
       IO_FLASH_INT, 1, &
       (/numFileBlks/), nodetype_str, len_trim(nodetype_str))
  call io_create_dataset(io_globalMe, fileID, libType, &
       IO_FLASH_INT, 1, &
       (/numFileBlks/), lrefine_str, len_trim(lrefine_str))
#ifdef FLASH_GRID_PARAMESH3OR4
  call io_create_dataset(io_globalMe, fileID, libType, &
       IO_FLASH_INT, 1, &
       (/numFileBlks/), which_child_str, len_trim(which_child_str))
  call io_create_dataset(io_globalMe, fileID, libType, &
       IO_FLASH_INT, 2, &
       (/numFileBlks,MFLAGS/), bflags_str, len_trim(bflags_str))
  call io_create_dataset(io_globalMe, fileID, libType, &
       IO_FLASH_INT, 5, &
       (/numFileBlks,1+(K3D*2),1+(K2D*2),1+(K1D*2),2/), &
       gr_gsurr_blks_str, len_trim(gr_gsurr_blks_str))
#endif
  call io_create_dataset(io_globalMe, fileID, libType, &
       IO_FLASH_INT, 2, &
       (/numFileBlks,nfaces+nchild+1/), gr_gid_str, len_trim(gr_gid_str))
  call io_create_dataset(io_globalMe, fileID, libType, &
       IO_FLASH_INT, 1, &
       (/numFileBlks/), procnumber_str, len_trim(procnumber_str))
  call io_create_dataset(io_globalMe, fileID, libType, &
       treeFP, 3, &
       (/numFileBlks,presentDims,2/), bnd_box_str, len_trim(bnd_box_str))
  call io_create_dataset(io_globalMe, fileID, libType, &
       treeFP, 2, &
       (/numFileBlks,presentDims/), coord_str, len_trim(coord_str))
  call io_create_dataset(io_globalMe, fileID, libType, &
       treeFP, 2, &
       (/numFileBlks,presentDims/), bsize_str, len_trim(bsize_str))

end subroutine io_createDatasets
