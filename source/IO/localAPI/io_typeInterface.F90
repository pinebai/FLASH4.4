!!****ih* source/IO/localAPI/io_typeInterface
!!
!! NAME
!!  io_typeInterface
!!
!! SYNOPSIS
!!  use io_typeInterface
!!
!! DESCRIPTION
!!
!! This is an interface module for Fortran subroutines that
!! are used in derived datatype I/O.
!!
!!***

#include "constants.h"
#include "Flash.h"
#include "io_flash.h"

module io_typeInterface
  interface
     subroutine io_xfer_mesh_data(fileID, fileFmt, fileType, &
     libType, xferType, localNumBlocks, globalBlockOffset)
       implicit none
       integer, intent(IN) :: fileID, fileFmt, fileType, libType, &
            xferType, localNumBlocks, globalBlockOffset
     end subroutine io_xfer_mesh_data
  end interface

  interface
     subroutine io_xfer_tree_data(tree_data, fileID, libType, xferType, &
          localNumBlocksIn, localOffsetIn, presentDims)
       use IO_data, ONLY : tree_data_t
       implicit none
       type(tree_data_t), intent(INOUT) :: tree_data
       integer, intent(IN) :: fileID, libType, xferType, &
            localNumBlocksIn, localOffsetIn, presentDims
     end subroutine io_xfer_tree_data
  end interface

  interface
     subroutine io_create_grid_header(myPE, fileID, fileFmt, fileType, &
          libType, dataFloatingPointType, metadataFloatingPointType)
       implicit none
       integer, intent(IN) :: myPE, fileID, fileFmt, fileType, &
            libType, dataFloatingPointType, metadataFloatingPointType
     end subroutine io_create_grid_header
  end interface

  interface
     subroutine io_getZeroBasedVarInfo(fileType, gridDataStruct, numGridVars, &
          numOutputGridVars, gridVarOffsets, gridVarLabels)
       implicit none
       integer, intent(IN) :: fileType, gridDataStruct
       integer, intent(OUT) :: numGridVars, numOutputGridVars 
       integer, dimension(MAX_MESH_VAR), intent(OUT) :: gridVarOffsets
       character (len=MAX_STRING_LENGTH), dimension(MAX_MESH_VAR), intent(OUT) :: &
            gridVarLabels
     end subroutine io_getZeroBasedVarInfo
  end interface

  interface
     subroutine io_getZeroBasedBlkSubarray(gridDataStruct, blockInnerSize, &
          blockOuterSize, blockInnerOffset)
       implicit none
       integer, intent(IN) :: gridDataStruct
       integer, dimension(MDIM), intent(OUT) :: blockInnerSize, &
            blockOuterSize, blockInnerOffset
     end subroutine io_getZeroBasedBlkSubarray
  end interface

  interface
     subroutine io_do_xfer(xferType, gridStruct, dataset, doXfer)
       implicit none
       integer, intent(IN) :: xferType, gridStruct
       character(len=*), intent(IN) :: dataset
       logical, intent(OUT) :: doXfer
     end subroutine io_do_xfer
  end interface
end module io_typeInterface
