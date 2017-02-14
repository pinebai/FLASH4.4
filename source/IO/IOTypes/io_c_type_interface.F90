!!****ih* source/IO/IOTypes/io_c_type_interface
!!
!! NAME
!!  io_c_type_interface
!!
!! SYNOPSIS
!!  use io_c_type_interface
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

module io_c_type_interface
  implicit none

#ifdef USE_IO_C_INTERFACE
  interface
     subroutine io_h5_read_file_format &
          (pMyPE, pFileID, pFileFormat) &
          bind(c)
       use iso_c_binding, only : c_int
       integer(c_int), intent(IN) :: pMyPE, pFileID
       integer(c_int), intent(OUT) :: pFileFormat
     end subroutine io_h5_read_file_format
  end interface

  interface
     subroutine io_xfer_mesh_dataset &
          (pMyPE, pFileID, pLibType, pXferType, pFileType, pFileFmt, &
          pGridDataStruct, pNumDataDims, pNumGridVars, pNonBlocking, &
          pPrePackData, diskOffset, memSize, memSubSize, &
          memOffset, memVarOffset, datasetName, pDsetNameLen, pData) &
          bind(c)
       use iso_c_binding, only : c_int, c_char, c_ptr
       integer(c_int), intent(IN) :: pMyPE, pFileID, pLibType, pXferType, &
            pFileType, pFileFmt, pGridDataStruct, pNumDataDims, pNumGridVars, &
            pNonBlocking, pPrePackData
       integer(c_int), dimension(*), intent(IN) :: diskOffset, &
            memSize, memSubSize, memOffset, memVarOffset
       character(kind=c_char), dimension(*), intent(IN) :: datasetName
       integer(c_int), intent(IN) :: pDsetNameLen
       type(c_ptr), value :: pData
     end subroutine io_xfer_mesh_dataset
  end interface

  interface
     subroutine io_init_grid_mpi_types &
          (pMyPE, pGridDataStruct, blockOuterSize, blockInnerSize, &
          blockInnerOffset, pNumGridVar, plotVarArr, pNumPlotVar) &
          bind(c)
       use iso_c_binding, only : c_int, c_char, c_ptr
       integer(c_int), intent(IN) :: pMyPE, pGridDataStruct
       integer(c_int), dimension(*), intent(IN) :: blockOuterSize, &
            blockInnerSize, blockInnerOffset
       integer(c_int), intent(IN) :: pNumGridVar
       integer(c_int), dimension(*), intent(IN) :: plotVarArr
       integer(c_int), intent(IN) :: pNumPlotVar
     end subroutine io_init_grid_mpi_types
  end interface
  
  interface
     subroutine io_free_grid_mpi_types() &
          bind(c)
     end subroutine io_free_grid_mpi_types
  end interface

  interface
     subroutine io_ncmpi_create_dimids &
          (pFileID, pDimVal, dimName, pDimNameLen) &
          bind(c)
       use iso_c_binding, only : c_int, c_char
       integer(c_int), intent(IN) :: pFileID, pDimVal
       character(kind=c_char), dimension(*), intent(IN) :: dimName
       integer(c_int), intent(IN) :: pDimNameLen
     end subroutine io_ncmpi_create_dimids
  end interface

  interface
     subroutine io_ncmpi_retrieve_dimids &
          (pMyPE, pFileID, pFileFmt, pGridStruct, dimIDs) &
          bind(c)
       use iso_c_binding, only : c_int
       integer(c_int), intent(IN) :: pMyPE, pFileID, pFileFmt, pGridStruct
       integer(c_int), dimension(*) :: dimIDs      
     end subroutine io_ncmpi_retrieve_dimids
  end interface

  interface
     subroutine io_ncmpi_nonblocking_init &
          (pMyPE) &
          bind(c)
       use iso_c_binding, only : c_int
       integer(c_int), intent(IN) :: pMyPE
     end subroutine io_ncmpi_nonblocking_init
  end interface

  interface
     subroutine io_ncmpi_nonblocking_complete_requests &
          (pFileID) &
          bind(c)
       use iso_c_binding, only : c_int
       integer(c_int), intent(IN) :: pFileID
     end subroutine io_ncmpi_nonblocking_complete_requests
  end interface

  interface
     subroutine io_ncmpi_nonblocking_finalize &
          () &
          bind(c)
     end subroutine io_ncmpi_nonblocking_finalize
  end interface
#endif

end module io_c_type_interface
