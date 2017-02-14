!!****if* source/IO/IOTypes/io_create_grid_header
!!
!! NAME
!!  io_create_grid_header
!!
!! SYNOPSIS
!!
!!  io_create_grid_header(integer(in) :: myPE,
!!                        integer(in) :: fileID,
!!                        integer(in) :: fileFmt,
!!                        integer(in) :: fileType,
!!                        integer(in) :: libType,
!!                        integer(in) :: dataFloatingPointType,
!!                        integer(in) :: metadataFloatingPointType)
!!                      
!!
!! DESCRIPTION
!!
!! This subroutine creates datasets (HDF5) and variables (pnetcdf)
!! for the mesh data structures e.g. UNK, FACEX.
!!
!! ARGUMENTS
!!
!! myPE: MPI process ID
!! fileID: The HDF5 or pnetcdf file identifier (used directly by the libraries)
!! fileFmt: The FLASH file layout.  The standard layout that tools
!!          understand is 9 (one mesh variable per dataset), but there is also
!!          support for an experimental file layout 10 (all mesh variables in
!!          the same dataset)
!! fileType: The FLASH file type (checkpoint file or plot file)
!! libType: The library we are using (HDF5 or pnetcdf)
!! dataFloatingPointType: Single precision or double precision for the
!!                        mesh data.
!! metadataFloatingPointType: Single precision or double precision for the
!!                            attributes of the mesh data.
!!
!!***

#include "constants.h"
#include "Flash.h"
#include "io_flash.h"

subroutine io_create_grid_header(myPE, fileID, fileFmt, fileType, &
     libType, dataFloatingPointType, metadataFloatingPointType)

#ifdef USE_IO_C_INTERFACE
  use iso_c_binding, ONLY : c_loc, c_char
  use io_c_interface, ONLY : io_create_dataset, io_attribute_create, &
       io_attribute_write, io_xfer_cont_slab
  use io_c_type_interface, ONLY : io_ncmpi_create_dimids, &
       io_ncmpi_retrieve_dimids
#endif
  use Grid_interface, ONLY : Grid_getGlobalIndexLimits
  use Driver_interface, ONLY : Driver_abortFlash  
  use io_typeInterface, ONLY : io_getZeroBasedVarInfo, &
       io_getZeroBasedBlkSubarray
  use IO_data, ONLY : io_splitNumBlks
  use io_typeData, ONLY : io_useLegacyLabels, io_legacyLabelLength

  implicit none
  integer, intent(IN) :: myPE, fileID, fileFmt, fileType, &
     libType, dataFloatingPointType, metadataFloatingPointType
#ifndef USE_IO_C_INTERFACE
#define c_loc(x) x
  integer, parameter :: c_char = KIND('A')
#endif

  integer, dimension(NUM_MESH_TYPES) :: gridDataStructs = &
       (/CENTER, SCRATCH, FACEX, FACEY, FACEZ/)

  !These are null terminated character arrays.  The c_loc memory address of
  !these arrays is interoperable with void *x.
  character(kind=c_char), dimension(MAX_MESH_VAR * MAX_STRING_LENGTH), target &
       :: packedGridLabels
  character(kind=c_char), dimension(MAX_STRING_LENGTH), target  :: packedVarStr

  !These are classic Fortran strings.  Interoperable with char *x or char x[].
  character(kind=c_char,len=MAX_STRING_LENGTH), dimension(MAX_MESH_VAR) :: &
       gridVarLabels
  character(kind=c_char,len=MAX_STRING_LENGTH) :: mesh_str, &
       num_var_str, var_str, strNum
  character(kind=c_char,len=7), parameter :: mesh_str_tbl(NUM_MESH_TYPES) = (/ &
       'unknown', &
       'facex  ', &
       'facey  ', &
       'facez  ', &
       'scratch' /)
  character(kind=c_char,len=19), parameter :: num_var_tbl(NUM_MESH_TYPES) = (/ &
       'dim_nUnkVar        ', &
       'dim_nFacexVar      ', &
       'dim_nFaceyVar      ', &
       'dim_nFacezVar      ', &
       'dim_nScratchGridVar' /)
  character(kind=c_char,len=9) :: nc_global_str = "NC_GLOBAL"
  character(kind=c_char,len=7) :: min_att_str = "minimum"
  character(kind=c_char,len=7) :: max_att_str = "maximum"


  integer, dimension(IO_MAX_DIMS) :: dimids
  integer, dimension(MAX_MESH_VAR) :: gridVarOffsets
  integer, dimension(MDIM) :: blockInnerSize, blockOuterSize, blockInnerOffset, &
       globalIndexLimits
  integer :: numGridVars, i, d, meshDims, labelLen, typeMatchedXfer, err
  integer, target :: numOutputGridVars

  interface
     subroutine pack_grid_labels(gridLabels, packedGridLabels, numMeshVar, strLen)
#ifdef USE_IO_C_INTERFACE
       use iso_c_binding, ONLY : c_char
#endif
       implicit none
#ifndef USE_IO_C_INTERFACE
#define c_loc(x) x
  integer, parameter :: c_char = KIND('A')
#endif
       character (kind=c_char,len=MAX_STRING_LENGTH), dimension(MAX_MESH_VAR), &
            intent(IN) :: gridLabels
       character (kind=c_char), dimension(MAX_MESH_VAR * MAX_STRING_LENGTH), &
            intent(OUT) :: packedGridLabels
       integer, intent(IN) :: numMeshVar, strLen
     end subroutine pack_grid_labels
  end interface
  interface
     subroutine create_null_char_array(a, b)
#ifdef USE_IO_C_INTERFACE
       use iso_c_binding, ONLY : c_char
#endif
       implicit none
#ifndef USE_IO_C_INTERFACE
#define c_loc(x) x
  integer, parameter :: c_char = KIND('A')
#endif
       character(kind=c_char, len=*), intent(IN) :: a
       character(kind=c_char), dimension(:), intent(OUT) :: b
     end subroutine create_null_char_array
  end interface

  !The io_xfer_cont_slab call transfers string data so casting from
  !double precision to single precision does not make sense.
  typeMatchedXfer = 0

  if (fileFmt <= 9) then
     meshDims = 4
  else if (fileFmt == 10) then
     meshDims = 5
  else
     call Driver_abortFlash("[io_create_grid_header]: Unknown file format")
  end if

  if (libType /= IO_FILE_PNETCDF .and. libType /= IO_FILE_HDF5) then
     call Driver_abortFlash("[io_create_grid_header]: Unknown library type")
  end if

  if (io_useLegacyLabels) then
     labelLen = io_legacyLabelLength
  else
     labelLen = MAX_STRING_LENGTH
  end if


  do d = lbound(gridDataStructs,1), ubound(gridDataStructs,1)

     call io_getZeroBasedVarInfo(fileType, gridDataStructs(d), numGridVars, &
          numOutputGridVars, gridVarOffsets, gridVarLabels)

     if (numOutputGridVars > 0) then

        !Total amount of global data to write: (Quite hacky at the moment).
        !---------------------------------------------------------------
        if (libType == IO_FILE_PNETCDF) then

#if defined(FLASH_IO_PNETCDF)
           num_var_str = trim(num_var_tbl(d))
           call io_ncmpi_create_dimids(fileID, numOutputGridVars, &
                num_var_str, len_trim(num_var_str))
           call io_ncmpi_retrieve_dimids(myPE, fileID, fileFmt, &
                gridDataStructs(d), dimids)
#endif

        else if (libType == IO_FILE_HDF5) then

#ifdef IO_FLASH_NOFBS_UG
           call Grid_getGlobalIndexLimits(globalIndexLimits)
           dimids(1) = 1
           dimids(2) = globalIndexLimits(3)
           dimids(3) = globalIndexLimits(2)
           dimids(4) = globalIndexLimits(1)
           dimids(5) = numOutputGridVars
           if (gridDataStructs(d) == SCRATCH) dimids(2:4) = dimids(2:4) + 1
           if (gridDataStructs(d) == FACEX) dimids(4) = dimids(4) + 1
           if (gridDataStructs(d) == FACEY) dimids(3) = dimids(3) + 1
           if (gridDataStructs(d) == FACEZ) dimids(2) = dimids(2) + 1
#else
           call io_getZeroBasedBlkSubarray(gridDataStructs(d), blockInnerSize, &
                blockOuterSize, blockInnerOffset)
           dimids(1) = io_splitNumBlks
           dimids(2) = blockInnerSize(3)
           dimids(3) = blockInnerSize(2)
           dimids(4) = blockInnerSize(1)
           dimids(5) = numOutputGridVars
#endif

        end if
        !---------------------------------------------------------------      
 

        !Create the datasets for each mesh data structure and the 
        !corresponding minimum and maximum attributes.  It is the same
        !for HDF5 and pnetcdf libraries.
        if (fileFmt <= 9) then
           do i = 1, numOutputGridVars
              var_str = trim(gridVarLabels(i))
              call io_create_dataset(myPE, fileID, libType, &
                   dataFloatingPointType, &
                   meshDims, dimids, var_str, len_trim(var_str))
              call io_attribute_create(myPE, fileID, libType, &
                   metadataFloatingPointType, &
                   1, (/1/), var_str, len_trim(var_str), &
                   min_att_str, len_trim(min_att_str))
              call io_attribute_create(myPE, fileID, libType, &
                   metadataFloatingPointType, &
                   1, (/1/), var_str, len_trim(var_str), &
                   max_att_str, len_trim(max_att_str))
           end do
        else
           mesh_str = trim(mesh_str_tbl(d))
           call io_create_dataset(myPE, fileID, libType, &
                dataFloatingPointType, &
                meshDims, dimids, mesh_str, len_trim(mesh_str))
           call io_attribute_create(myPE, fileID, libType, &
                metadataFloatingPointType, &
                1, (/numOutputGridVars/), mesh_str, &
                len_trim(mesh_str), min_att_str, len_trim(min_att_str))
           call io_attribute_create(myPE, fileID, libType, &
                metadataFloatingPointType, &
                1, (/numOutputGridVars/), mesh_str, &
                len_trim(mesh_str), max_att_str, len_trim(max_att_str))
        end if


        !Create datasets for the labels describing the mesh data structures.
        if (libType == IO_FILE_PNETCDF) then
           mesh_str = trim(mesh_str_tbl(d))//"_names"

           call io_attribute_create(myPE, fileID, libType, &
                IO_FLASH_INT, 1, (/1/), nc_global_str, len_trim(nc_global_str), &
                mesh_str, len_trim(mesh_str))
           call io_attribute_write(myPE, fileID, libType, &
                IO_FLASH_INT, nc_global_str, len_trim(nc_global_str), &
                mesh_str, len_trim(mesh_str), &
                c_loc(numOutputGridVars))
           do i = 1, numOutputGridVars
              Write(strNum,'(i10)') i-1
              strNum = adjustl(strNum)
              mesh_str = trim(mesh_str_tbl(d))//"_"//trim(strNum)
              call io_attribute_create(myPE, fileID, libType, &
                   IO_FLASH_CHAR, 1, (/labelLen/), nc_global_str, &
                   len_trim(nc_global_str), mesh_str, &
                   len_trim(mesh_str))
              call create_null_char_array(gridVarLabels(i), packedVarStr)
              call io_attribute_write(myPE, fileID, libType, &
                   IO_FLASH_CHAR, nc_global_str, len_trim(nc_global_str), &
                   mesh_str, len_trim(mesh_str), &
                   c_loc(packedVarStr(1)))
           end do
        else if (libType == IO_FILE_HDF5) then
           mesh_str = trim(mesh_str_tbl(d))//" names"

           call pack_grid_labels(gridVarLabels, packedGridLabels, &
                numOutputGridVars, labelLen)
           call io_create_dataset(myPE, fileID, libType, &
                IO_FLASH_STRING, &
                3, (/numOutputGridVars,1,labelLen/), mesh_str, &
                len_trim(mesh_str))
           call io_xfer_cont_slab(myPE, fileID, IO_FILE_HDF5, &
                IO_WRITE_XFER_MASTER_PE, typeMatchedXfer, mesh_str, &
                len_trim(mesh_str), IO_FLASH_STRING, &
                (/numOutputGridVars,1,labelLen/), (/0,0,0/), &
                (/numOutputGridVars,1,labelLen/), (/0,0,0/), &
                (/numOutputGridVars,1,labelLen/), 3, c_loc(packedGridLabels(1)), &
                err)
        end if
     end if
  end do
end subroutine io_create_grid_header


subroutine pack_grid_labels(gridLabels, packedGridLabels, numMeshVar, strLen)

#ifdef USE_IO_C_INTERFACE
  use iso_c_binding, ONLY : c_char
#endif

  implicit none
#ifndef USE_IO_C_INTERFACE
#define c_loc(x) x
  integer, parameter :: c_char = KIND('A')
#endif
  character (kind=c_char,len=MAX_STRING_LENGTH), dimension(MAX_MESH_VAR), &
       intent(IN) :: gridLabels
  character (kind=c_char), dimension(MAX_MESH_VAR * MAX_STRING_LENGTH), &
       intent(OUT) :: packedGridLabels
  integer, intent(IN) :: numMeshVar, strLen
  integer :: i, j, k
  character (kind=c_char), dimension(MAX_STRING_LENGTH) :: tmpLabel
  interface
     subroutine create_null_char_array(a, b)
#ifdef USE_IO_C_INTERFACE
       use iso_c_binding, ONLY : c_char
#endif
       implicit none
#ifndef USE_IO_C_INTERFACE
#define c_loc(x) x
  integer, parameter :: c_char = KIND('A')
#endif
       character(kind=c_char, len=*), intent(IN) :: a
       character(kind=c_char), dimension(:), intent(OUT) :: b
     end subroutine create_null_char_array
  end interface

  !Initialize all packed data to keep valgrind happy.
  do i = 1, MAX_MESH_VAR * MAX_STRING_LENGTH
     packedGridLabels(i) = CHAR(0)
  end do

  k = 1
  do i = 1, numMeshVar
     call create_null_char_array(gridLabels(i), tmpLabel)
     do j = 1, strLen
        !'dens', 'pres', 'temp' labels will be stored as 'densprestemp'.
        packedGridLabels(k) = tmpLabel(j)
        k = k + 1
     end do
  end do
end subroutine pack_grid_labels


subroutine create_null_char_array(a, b)

  use Driver_interface, ONLY : Driver_abortFlash
#ifdef USE_IO_C_INTERFACE
  use iso_c_binding, ONLY : c_char
#endif

  implicit none
#ifndef USE_IO_C_INTERFACE
#define c_loc(x) x
  integer, parameter :: c_char = KIND('A')
#endif
  character(kind=c_char, len=*), intent(IN) :: a
  character(kind=c_char), dimension(:), intent(OUT) :: b
  integer :: i, aLen, bLen, copyLen

  !We assume a simple unit based character array.
  if (lbound(b,1) /= 1) then
     call Driver_abortFlash("[create_null_char_array]: Array not unit based")
  end if

  !aLen contains the number of characters excluding trailing blanks in a.
  !bLen contain the number of characters in the memory space of b.
  aLen = len_trim(a)
  bLen = size(b, 1)

  !Copy the trimmed string.
  copyLen = min(aLen, bLen)

  do i = 1, copyLen
     b(i) = a(i:i)
  end do

  if (bLen > copyLen) then
     !Initialize all remaining space so that valgrind does not complain
     !when uninitialized memory is copied around in e.g. HDF5 library.
     do i = copyLen + 1, bLen
        b(i) = CHAR(0)
     end do
  else
     !Null terminate the string even if it means losing a valid character.
     b(bLen) = CHAR(0)
  end if
end subroutine create_null_char_array
