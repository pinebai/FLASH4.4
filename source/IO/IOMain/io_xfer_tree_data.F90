!!****if* source/IO/IOMain/io_xfer_tree_data
!!
!! NAME
!!  io_xfer_tree_data
!!
!! SYNOPSIS
!!
!!  call io_xfer_tree_data(type(tree_data_t)(INOUT) :: tree_data,
!!                         integer(in)              :: fileID,
!!                         integer(in)              :: libType,
!!                         integer(in)              :: xferType,
!!                         integer(in)              :: localNumBlocks,
!!                         integer(in)              :: localOffset,
!!                         integer(in)              :: presentDims)
!!                      
!!
!! DESCRIPTION
!!
!! This subroutine transfers tree data such as block refinement level, block
!! bounding box and block size from
!!   1. memory to file for a write transfer, e.g. write checkpoint.
!!   2. file to memory for a read transfer, e.g. read checkpoint.
!!
!! The tree data is either read from or written to a data structure named
!! tree_data.  In the case of a read transfer it is assumed that the fields
!! of tree_data have been pre-allocated - we assert that the pre-allocated
!! storage is large enough in the C data transfer function io_xfer_cont_slab.
!!
!! Note that tree data does not need to be read from file in UG applications.
!!
!!
!! ARGUMENTS
!!
!! tree_data: The actual tree data.
!! fileID: The HDF5 or pnetcdf file identifier (used directly by the libraries).
!! libType: The library we are using (HDF5 or pnetcdf).
!! xferType: The direction of data transfer: memory to file or file to memory.
!! localNumBlocks: The number of blocks on myPE being transferred to/from file.
!! localOffset: The read/write block offset in file.
!! presentDims: The number of dimensions in the coord, bsize and
!!              bndbox datasets.  Ensures backwards compatibility with
!!              FLASH3 beta and earlier.
!!***


#include "constants.h"
#include "Flash.h"
#include "io_flash.h"

subroutine io_xfer_tree_data(tree_data, fileID, &
     libType, xferType, localNumBlocks, localOffset, presentDims)

#ifdef USE_IO_C_INTERFACE
  use iso_c_binding, ONLY : c_loc, c_char
  use io_c_interface, ONLY : io_xfer_cont_slab
#endif

#ifdef FLASH_GRID_PARAMESH
  use tree, ONLY : nfaces, nchild, lrefine_max
# ifdef FLASH_GRID_PARAMESH3OR4
  use tree, ONLY : MFLAGS
  use Grid_data, ONLY : gr_is_gsurr_blks_initialized
# endif
#endif

  use Driver_interface, ONLY : Driver_abortFlash
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use IO_data, ONLY : io_type_matched_xfer, io_globalMe, tree_data_t

  implicit none
  type(tree_data_t), intent(INOUT) :: tree_data
  integer, intent(IN) :: fileID, libType, xferType, localNumBlocks, &
       localOffset, presentDims

#ifndef USE_IO_C_INTERFACE
#define c_loc(x) x
  integer, parameter :: c_char = KIND('A')
#endif

  !All strings must be null-terminated.
  character(kind=c_char,len=MAX_STRING_LENGTH) :: bflags_str, gr_gid_str, &
       coord_str, nodetype_str, lrefine_str, which_child_str, &
       procnumber_str, bnd_box_str, bsize_str, tree_str, gr_gsurr_blks_str

#ifndef FLASH_GRID_PARAMESH
  integer, parameter :: nfaces = 2*NDIM, nchild = 2**NDIM, lrefine_max = 1
#endif
  integer :: typeMatchedXfer, err
  logical :: doDataXfer
#ifdef FLASH_GRID_PARAMESH4DEV_SURR_BLKS_OPTION
  logical, parameter :: do_gsurr_blks_read = .true.
#else
  logical, parameter :: do_gsurr_blks_read = .false.
#endif

  if (xferType == IO_WRITE_XFER .or. xferType == IO_WRITE_XFER_MASTER_PE) then
     tree_str = "write tree data"
  else if (xferType == IO_READ_XFER .or. xferType == IO_READ_XFER_MASTER_PE) then
     tree_str = "read tree data"
  else
     call Driver_abortFlash("Invalid transfer type")
  end if


  doDataXfer = .true.
#ifdef FLASH_GRID_UG
  if (xferType == IO_READ_XFER .or. xferType == IO_READ_XFER_MASTER_PE) then
     doDataXfer = .false.
  end if
#endif


  if (doDataXfer) then
     IO_TIMERS_START(tree_str)

     !This option ensures that file reads and writes involve floating point
     !data of the same type in both memory and file.  It is essential to
     !type match in large parallel runs that use collective I/O with HDF5.
     typeMatchedXfer = 0
     if (libType == IO_FILE_HDF5) then
        if (io_type_matched_xfer) then
           typeMatchedXfer = 1
        end if
     end if

     bflags_str = "bflags"
     gr_gid_str = "gid"
     gr_gsurr_blks_str = "gsurr_blks"
     coord_str = "coordinates"

     if (libType == IO_FILE_HDF5) then
        nodetype_str = "node type"
        lrefine_str = "refine level"
        which_child_str = "which child"
        procnumber_str = "processor number"
        bnd_box_str = "bounding box"
        bsize_str = "block size"
     else if (libType == IO_FILE_PNETCDF) then
        nodetype_str = "nodetype"
        lrefine_str = "lrefine"
        which_child_str = "which_child"
        procnumber_str = "processor_number"
        bnd_box_str = "bndbox"
        bsize_str = "blocksize"
     else
        call Driver_abortFlash("[io_xfer_tree_data]: Unrecognized file type")
     end if



     IO_TIMERS_START(nodetype_str)
     if (.not.associated(tree_data % nodetype)) then
        call Driver_abortFlash("[io_xfer_tree_data]: nodetype not associated")
     end if
     !Paramesh memory size: (/maxblocks_tr/)
     call io_xfer_cont_slab(io_globalMe, fileID, libType, &
          xferType, typeMatchedXfer, nodetype_str, len_trim(nodetype_str), &
          IO_FLASH_INT, (/size(tree_data % nodetype,1)/), (/0/), &
          (/localNumBlocks/), (/localOffset/), &
          (/localNumBlocks/), 1, c_loc(tree_data % nodetype(1)), err)
     IO_CHECK_XFER(err, nodetype_str)
     IO_TIMERS_STOP(nodetype_str)


     IO_TIMERS_START(lrefine_str)
     if (.not.associated(tree_data % lrefine)) then
        call Driver_abortFlash("[io_xfer_tree_data]: lrefine not associated")
     end if
     !Paramesh memory size: (/maxblocks_tr/)
     call io_xfer_cont_slab(io_globalMe, fileID, libType, &
          xferType, typeMatchedXfer, lrefine_str, len_trim(lrefine_str), &
          IO_FLASH_INT, (/size(tree_data % lrefine,1)/), (/0/), &
          (/localNumBlocks/), (/localOffset/), &
          (/localNumBlocks/), 1, c_loc(tree_data % lrefine(1)), err)
     IO_CHECK_XFER(err, lrefine_str)
     IO_TIMERS_STOP(lrefine_str)
     if (xferType == IO_READ_XFER .or. xferType == IO_READ_XFER_MASTER_PE) then
        !Abort if the checkpoint file contains blocks at refinement
        !level greater than lrefine_max.  This prevents bad things
        !happening later on.
        if (any(tree_data % lrefine(1:localNumBlocks) > lrefine_max)) then
           write(6,'(2(a,i4))') &
                " Max block refinement level in the checkpoint file ", &
                maxval(tree_data % lrefine(1:localNumBlocks)), &
                ", lrefine_max ", lrefine_max
           call Driver_abortFlash("At least 1 block in the checkpoint "// &
                "file has a refinement level > lrefine_max. "// &
                "Increase lrefine_max in your flash.par!")
        end if
     end if


#ifdef FLASH_GRID_PARAMESH3OR4
     !data struct new to PM3.  which_child definitely
     !needed for restart. Adding bflags for completeness.
     IO_TIMERS_START(which_child_str)
     if (.not.associated(tree_data % which_child)) then
        call Driver_abortFlash("[io_xfer_tree_data]: which_child not associated")
     end if
     !Paramesh memory size: (/maxblocks_tr/)
     call io_xfer_cont_slab(io_globalMe, fileID, libType, &
          xferType, typeMatchedXfer, which_child_str, &
          len_trim(which_child_str), &
          IO_FLASH_INT, (/size(tree_data % which_child,1)/), (/0/), &
          (/localNumBlocks/), (/localOffset/), &
          (/localNumBlocks/), 1, c_loc(tree_data % which_child(1)), err)
     IO_CHECK_XFER(err, which_child_str)
     IO_TIMERS_STOP(which_child_str)


     if (libType == IO_FILE_HDF5) then
        IO_TIMERS_START(bflags_str)
        if (.not.associated(tree_data % bflags)) then
           call Driver_abortFlash("[io_xfer_tree_data]: bflags not associated")
        end if
        !Paramesh memory size: (/maxblocks_tr,MFLAGS/)
        call io_xfer_cont_slab(io_globalMe, fileID, libType, &
             xferType, typeMatchedXfer, bflags_str, len_trim(bflags_str), &
             IO_FLASH_INT, (/size(tree_data % bflags,2), &
             size(tree_data % bflags,1)/), (/0,0/), &
             (/localNumBlocks,MFLAGS/), (/localOffset,0/), &
             (/localNumBlocks,MFLAGS/), 2, c_loc(tree_data % bflags(1,1)), err)
        IO_CHECK_XFER(err, bflags_str)
        IO_TIMERS_STOP(bflags_str)
     end if


     if ((xferType == IO_WRITE_XFER .or. xferType == IO_WRITE_XFER_MASTER_PE) .or. &
          ((xferType == IO_READ_XFER .or. xferType == IO_READ_XFER_MASTER_PE) .and. &
          do_gsurr_blks_read)) then
        IO_TIMERS_START(gr_gsurr_blks_str)
        if (.not.associated(tree_data % gsurr_blks)) then
           call Driver_abortFlash("[io_xfer_tree_data]: gsurr_blks not associated")
        end if
        call io_xfer_cont_slab(io_globalMe, fileID, libType, &
             xferType, typeMatchedXfer, &
             gr_gsurr_blks_str, len_trim(gr_gsurr_blks_str), &
             IO_FLASH_INT, (/size(tree_data % gsurr_blks,5), &
             size(tree_data % gsurr_blks,4),size(tree_data % gsurr_blks,3), &
             size(tree_data % gsurr_blks,2),size(tree_data % gsurr_blks,1)/), &
             (/0,0,0,0,0/), &
             (/localNumBlocks,1+(K3D*2),1+(K2D*2),1+(K1D*2),2/), &
             (/localOffset,0,0,0,0/), &
             (/localNumBlocks,1+(K3D*2),1+(K2D*2),1+(K1D*2),2/), &
             5, c_loc(tree_data % gsurr_blks(1,1,1,1,1)), err)
        if (xferType == IO_READ_XFER .or. xferType == IO_READ_XFER_MASTER_PE) then
           gr_is_gsurr_blks_initialized = (err == 0)
        else
           IO_CHECK_XFER(err, gr_gsurr_blks_str)
        end if
        IO_TIMERS_STOP(gr_gsurr_blks_str)
     else
        gr_is_gsurr_blks_initialized = .false.
     end if
#endif


     IO_TIMERS_START(gr_gid_str)
     if (.not.associated(tree_data % gid)) then
        call Driver_abortFlash("[io_xfer_tree_data]: gid not associated")
     end if
     !Paramesh memory size: (/MAXBLOCKS,nfaces+nchild+1/)
     call io_xfer_cont_slab(io_globalMe, fileID, libType, &
          xferType, typeMatchedXfer, gr_gid_str, len_trim(gr_gid_str), &
          IO_FLASH_INT, (/size(tree_data % gid,2), &
          size(tree_data % gid,1)/), (/0,0/), &
          (/localNumBlocks,nfaces+nchild+1/), (/localOffset,0/), &
          (/localNumBlocks,nfaces+nchild+1/), 2, c_loc(tree_data % gid(1,1)), &
          err)
     IO_CHECK_XFER(err, gr_gid_str)
     IO_TIMERS_STOP(gr_gid_str)


     if (xferType == IO_WRITE_XFER .or. xferType == IO_WRITE_XFER_MASTER_PE) then
        IO_TIMERS_START(procnumber_str)
        if (.not.associated(tree_data % procnumber)) then
           call Driver_abortFlash("[io_xfer_tree_data]: procnumber not associated")
        end if
        !Paramesh memory size: (/localNumBlocks/)
        call io_xfer_cont_slab(io_globalMe, fileID, libType, &
             xferType, typeMatchedXfer, procnumber_str, &
             len_trim(procnumber_str), &
             IO_FLASH_INT,(/size(tree_data % procnumber,1)/), (/0/), &
             (/localNumBlocks/), (/localOffset/), &
             (/localNumBlocks/), 1, c_loc(tree_data % procnumber(1)), err)
        IO_CHECK_XFER(err, procnumber_str)
        IO_TIMERS_STOP(procnumber_str)
     end if


     IO_TIMERS_START(bnd_box_str)
     if (.not.associated(tree_data % bnd_box)) then
        call Driver_abortFlash("[io_xfer_tree_data]: bnd_box not associated")
     end if
     !Paramesh memory size: (/maxblocks_tr,MDIM,2/)
     call io_xfer_cont_slab(io_globalMe, fileID, libType, &
          xferType, typeMatchedXfer, bnd_box_str, len_trim(bnd_box_str), &
          IO_FLASH_DOUBLE, (/size(tree_data % bnd_box,3),&
          size(tree_data % bnd_box,2),size(tree_data % bnd_box,1)/), (/0,0,0/), &
          (/localNumBlocks,presentDims,2/), (/localOffset,0,0/), &
          (/localNumBlocks,presentDims,2/), 3, &
          c_loc(tree_data % bnd_box(1,1,1)), err)
     IO_CHECK_XFER(err, bnd_box_str)
     IO_TIMERS_STOP(bnd_box_str)


     !Block center coordinates (not to be confused with the cell coordinates)
     IO_TIMERS_START(coord_str)
     if (.not.associated(tree_data % coord)) then
        call Driver_abortFlash("[io_xfer_tree_data]: coord not associated")
     end if
     !Paramesh memory size: (/maxblocks_tr,MDIM/)
     call io_xfer_cont_slab(io_globalMe, fileID, libType, &
          xferType, typeMatchedXfer, coord_str, len_trim(coord_str), &
          IO_FLASH_DOUBLE, (/size(tree_data % coord,2), &
          size(tree_data % coord,1)/), (/0,0/), &
          (/localNumBlocks,presentDims/), (/localOffset,0/), &
          (/localNumBlocks,presentDims/), 2, &
          c_loc(tree_data % coord(1,1)), err)
     IO_CHECK_XFER(err, coord_str)
     IO_TIMERS_STOP(coord_str)


     IO_TIMERS_START(bsize_str)
     if (.not.associated(tree_data % bsize)) then
        call Driver_abortFlash("[io_xfer_tree_data]: bsize not associated")
     end if
     !Paramesh memory size: (/maxblocks_tr,MDIM/)
     call io_xfer_cont_slab(io_globalMe, fileID, libType, &
          xferType, typeMatchedXfer, bsize_str, len_trim(bsize_str), &
          IO_FLASH_DOUBLE, (/size(tree_data % bsize,2), &
          size(tree_data % bsize,1)/), (/0,0/), &
          (/localNumBlocks,presentDims/), (/localOffset,0/), &
          (/localNumBlocks,presentDims/), 2, &
          c_loc(tree_data % bsize(1,1)), err)
     IO_CHECK_XFER(err, bsize_str)
     IO_TIMERS_STOP(bsize_str)


     IO_TIMERS_STOP(tree_str)
  end if
end subroutine io_xfer_tree_data
