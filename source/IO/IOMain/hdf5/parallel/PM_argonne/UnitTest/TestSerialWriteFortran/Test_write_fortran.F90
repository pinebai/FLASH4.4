!!****if* source/IO/IOMain/hdf5/parallel/PM_argonne/UnitTest/TestSerialWriteFortran/Test_write_fortran
!!
!! NAME
!!
!! Test_write_fortran
!!
!!
!! SYNOPSIS
!!
!! Test_write_fortran()
!!
!! DESCRIPTION
!!
!! This is a standalone program that tests I/O independently of FLASH
!!
!! ARGUMENTS
!!
!! NOTES
!!
!! We assume code is compiled with an option that promotes 
!! 4 bytes reals to 8 byte reals.
!!
!!
!! We want to write the internal region of dataArray to file.
!! This is: dataArray(1, 2:3, 2:4, 1, 1) as NXB=2, NYB=3, NGUARD=1.
!!
!!
!! In Fortran the left most dimension in a multi-dimensional array is 
!! contiguous.  So for an array(i,j), I like to draw the contiguous axis (i)
!! on the column because Fortran is column major.
!!
!!  i |
!!    |
!!    -----
!!      j
!!
!!              MEMORY-SPACE                  DISK-SPACE
!!
!!      i=1 | 01  05  09  13  17     i=1 | **  **  **  **  **
!!      i=2 | 02  06  10  14  18     i=2 | **  06  10  14  **
!!      i=3 | 03  07  11  15  19     i=3 | **  07  11  15  **
!!      i=4 | 04  08  12  16  20     i=4 | **  **  **  **  **
!!          --------------------         --------------------
!!           j=1 j=2 j=3 j=4 j=5          j=1 j=2 j=3 j=4 j=5
!!
!!
!! In this unit-test the HDF5 file should contain 6,7,10,11,14,15.
!!
!!***

Program Test_write_fortran

#include "constants.h"
#include "io_flash.h"

!We will have this many variables:
#define NUNK_VARS 1
#define VAR_LEN 4

implicit none


!Metadata definitions / declarations:
!-------------------------------------------------------------------------------
real, parameter, dimension(NUNK_VARS) :: minVar = (/ 10.1 /)
real, parameter, dimension(NUNK_VARS) :: maxVar = (/ 22.1 /)
integer, parameter :: numUnkVar = NUNK_VARS
character(len=*), parameter :: attMinName = "minimum"
character(len=*), parameter :: attMaxName = "maximum"

!Data definitions / declarations:
!-------------------------------------------------------------------------------
integer, parameter :: NXB = 2, NYB = 3, NGUARD = 1
real, dimension(1, NXB+2*NGUARD, NYB+2*NGUARD, 1, 1) :: dataArray

!These arrays are how I handle both SCRATCH and UNK layouts which have 
!NXB, NYB, NZB in different array slots.
integer, parameter, dimension(IO_MESH_DIMS) :: &
     numInternalCells = (/ 1, NXB, NYB, 1, 1 /)
integer, parameter, dimension(IO_MESH_DIMS) :: &
     numGuardCells = (/ 0, NGUARD, NGUARD, 0, 0 /)
integer, dimension(IO_MESH_DIMS) :: numInternalCells_C, numGuardCells_C, &
     globalSize, globalOffset, localSize, localSubSize, localOffset
integer, parameter :: numDataDims = IO_MESH_DIMS

!Other definitions / declarations:
!-------------------------------------------------------------------------------
integer :: i, j, p, count, fileID, err
integer, parameter :: myPE = 0
character (len=*), parameter :: datasetName = "unknown"
integer :: typeMatchedXfer = 1

external init_hdf5_c, finalise_hdf5_c, io_create_dataset, io_xfer_cont_slab, &
     io_attribute_create, io_attribute_write


!Prepare the data.
!-------------------------------------------------------------------------------
count = 1
do j = 1, NYB+2*NGUARD
   do i = 1, NXB+2*NGUARD
      dataArray(1,i,j,1,1) = real(count)
      count = count + 1
   end do
end do
print *, "Single process test: HDF5 file should contain 6,7,10,11,14,15"


!Reverse the array of indicies into a C-order.
do p = 1, numDataDims
   numInternalCells_C(p) = numInternalCells((numDataDims+1) - p)
   numGuardCells_C(p) = numGuardCells((numDataDims+1) - p)
end do
globalSize = numInternalCells_C
globalOffset = 0
localSubSize = numInternalCells_C
localSize = numInternalCells_C + (2 * numGuardCells_C)
localOffset = numGuardCells_C


!Initialise HDF5.  We get back fileID.
call init_hdf5_c(fileID)


call io_create_dataset(myPE, &
     fileID, &
     IO_FILE_HDF5, &
     IO_FLASH_DOUBLE, &
     numDataDims, &
     globalSize, &
     datasetName, &
     len_trim(datasetName))

call io_attribute_create(myPE, &
     fileID, &
     IO_FILE_HDF5, &
     IO_FLASH_DOUBLE, &
     numUnkVar, &
     numUnkVar, &
     datasetName, &
     len_trim(datasetName), &
     attMinName, &
     len_trim(attMinName))

call io_attribute_create(myPE, &
     fileID, &
     IO_FILE_HDF5, &
     IO_FLASH_DOUBLE, &
     numUnkVar, &
     numUnkVar, &
     datasetName, &
     len_trim(datasetName), &
     attMaxName, &
     len_trim(attMaxName))

call io_xfer_cont_slab(myPE, &
     fileID, &
     IO_FILE_HDF5, &
     IO_WRITE_XFER, &
     typeMatchedXfer, &
     datasetName, &
     len_trim(datasetName), &
     IO_FLASH_DOUBLE, &
     localSize, &
     localOffset, &
     localSubSize, &
     globalOffset, &
     localSubSize, &
     numDataDims, &
     dataArray, &
     err)
if (err /= 0) then
   print *, "I/O error"
   stop
end if

call io_attribute_write(myPE, &
     fileID, &
     IO_FILE_HDF5, &
     IO_FLASH_DOUBLE, &
     datasetName, &
     len_trim(datasetName), &
     attMinName, &
     len_trim(attMinName), &
     minVar)

call io_attribute_write(myPE, &
     fileID, &
     IO_FILE_HDF5, &
     IO_FLASH_DOUBLE, &
     datasetName, &
     len_trim(datasetName), &
     attMaxName, &
     len_trim(attMaxName), &
     maxVar)

call finalise_hdf5_c(fileID)

End Program Test_write_fortran
