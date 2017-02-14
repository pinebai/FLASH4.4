!!****ih* source/IO/IOMain/io_c_interface
!!
!! NAME
!!  io_c_interface
!!
!! SYNOPSIS
!!  use io_c_interface
!!
!! DESCRIPTION
!!
!! This is an interface module for C I/O functions.  It is not
!! essential because our Fortran subroutines only ever pass
!! primitive types, however, it offers type checking and
!! guaranteed interoperability of strings.
!! 
!! It is conditionally compiled because old compilers do not 
!! support Fortran 2003 features.
!!
!! Interoperability of void pointers:
!! type(c_ptr), value -> void *
!! type(c_ptr) -> void **
!!
!!***

#include "constants.h"
#include "Flash.h"

module io_c_interface
  implicit none

#ifdef USE_IO_C_INTERFACE
  interface
     subroutine io_create_dataset &
          (pMyPE, pFileID, pLibType, pDiskType, pDims, dimIDs, datasetName, &
          pDsetNameLen) &
          bind(c)
       use iso_c_binding, only : c_int, c_char
       integer(c_int), intent(IN) :: pMyPE, pFileID, pLibType, pDiskType, pDims
       integer(c_int), dimension(*), intent(IN) :: dimIDs
       character(kind=c_char), dimension(*), intent(IN) :: datasetName
       integer(c_int), intent(IN) :: pDsetNameLen
     end subroutine io_create_dataset
  end interface

  interface
     subroutine io_attribute_create &
          (pMyPE, pFileID, pLibType, pDiskType, pDims, datasetSize, &
          datasetName, pDsetNameLen, attDatasetName, pAttNameLen) &
          bind(c)
       use iso_c_binding, only : c_int, c_char
       integer(c_int), intent(IN) :: pMyPE, pFileID, pLibType, pDiskType, pDims
       integer(c_int), dimension(*), intent(IN) :: datasetSize
       character(kind=c_char), dimension(*), intent(IN) :: datasetName
       integer(c_int), intent(IN) :: pDsetNameLen
       character(kind=c_char), dimension(*), intent(IN) :: attDatasetName
       integer(c_int), intent(IN) :: pAttNameLen
     end subroutine io_attribute_create
  end interface

  interface
     subroutine io_attribute_write &
          (pMyPE, pFileID, pLibType, pMemType, datasetName, &
          pDsetNameLen, attDatasetName, pAttNameLen, pData) &
          bind(c)
       use iso_c_binding, only : c_int, c_char, c_ptr
       integer(c_int), intent(IN) :: pMyPE, pFileID, pLibType, pMemType
       character(kind=c_char), dimension(*), intent(IN) :: datasetName
       integer(c_int), intent(IN) :: pDsetNameLen
       character(kind=c_char), dimension(*), intent(IN) :: attDatasetName
       integer(c_int), intent(IN) :: pAttNameLen
       type(c_ptr), value, intent(IN) :: pData
     end subroutine io_attribute_write
  end interface

  interface
     subroutine io_xfer_cont_slab &
          (pMyPE, pFileID, pFileType, pXferType, pTypeMatchedXfer, &
          datasetName, pNameLength, &
          pMemType, memSize, memStart, memCount, diskStart, diskCount, &
          pDims, pData, pErr) &
          bind(c)
       use iso_c_binding, only : c_int, c_char, c_ptr
       integer(c_int), intent(IN) :: pMyPE, pFileID, pFileType, pXferType, &
            pTypeMatchedXfer
       character(kind=c_char), dimension(*), intent(IN) :: datasetName
       integer(c_int), intent(IN) :: pNameLength, pMemType
       integer(c_int), dimension(*), intent(IN) :: memSize, memStart, &
            memCount, diskStart, diskCount
       integer(c_int), intent(IN) :: pDims
       type(c_ptr), value :: pData
       integer(c_int), intent(OUT) :: pErr
     end subroutine io_xfer_cont_slab
  end interface

  interface
     subroutine io_ncmpi_name_to_id &
          (pFileID,name,pLen,pVarID) &
          bind(c)
       use iso_c_binding, only : c_int, c_char
       integer(c_int), intent(IN) :: pFileID
       character(kind=c_char), dimension(*), intent(IN) :: name
       integer(c_int), intent(IN) :: pLen
       integer(c_int), intent(OUT) :: pVarID
     end subroutine io_ncmpi_name_to_id
  end interface

  interface
     subroutine io_ncmpi_define_mode_redef &
          (pFileID) &
          bind(c)
       use iso_c_binding, only : c_int
       integer(c_int), intent(IN) :: pFileID
     end subroutine io_ncmpi_define_mode_redef
  end interface

  interface
     subroutine io_ncmpi_define_mode_enddef &
          (pFileID) &
          bind(c)
       use iso_c_binding, only : c_int
       integer(c_int), intent(IN) :: pFileID
     end subroutine io_ncmpi_define_mode_enddef
  end interface

  interface
     subroutine io_ncmpi_read_file_format &
          (pMyPE,pFileID,pFileFormat) &
          bind(c)
       use iso_c_binding, only : c_int
       integer(c_int), intent(IN) :: pMyPE,pFileID
       integer(c_int), intent(OUT) :: pFileFormat
     end subroutine io_ncmpi_read_file_format
  end interface
#endif

end module io_c_interface
