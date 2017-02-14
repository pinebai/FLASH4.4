!!****if* source/IO/IOTypes/io_typeInit
!!
!! NAME
!!  io_typeInit
!!
!! SYNOPSIS
!!
!!  io_typeInit(integer(in) :: myPE)
!!              integer(in) :: numProcs)
!!
!! DESCRIPTION
!!
!!  Perform initialization for IO Types subunit.
!!
!! ARGUMENTS
!!
!!  myPE: MPI process identifier
!!  numProcs:  Total number of processes.
!!
!!***

#include "constants.h"
#include "Flash.h"
#include "io_flash.h"

subroutine io_typeInit(myPE, numProcs)
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_abortFlash
  use io_typeInterface, ONLY : io_getZeroBasedVarInfo, &
       io_getZeroBasedBlkSubarray
#ifdef USE_IO_C_INTERFACE
  use io_c_type_interface, ONLY : io_init_grid_mpi_types
#ifdef FLASH_IO_PNETCDF
  use io_c_type_interface, ONLY : io_ncmpi_nonblocking_init
#endif
#endif
  use Grid_interface, ONLY : Grid_getNumVars
  use IO_data, ONLY : io_unklabelsGlobal, io_fileFormatVersion
  use io_typeData

  implicit none
  include "Flash_mpi.h"

  integer, intent(IN) :: myPE, numProcs
  integer, dimension(5) :: allGridDataStruct = &
       (/CENTER, SCRATCH, FACEX, FACEY, FACEZ/)
  !There can't be more plot variables than checkpoint variables.
  integer, dimension(MAX_MESH_VAR) :: plotVarOffsets
  integer, dimension(MDIM) :: blockInnerSize, blockOuterSize, blockInnerOffset
  integer :: i, numGridVars, numPlotVars, gridDataStruct, meshCopyCount, &
       extentVarDim
  character (len=MAX_STRING_LENGTH), dimension(MAX_MESH_VAR) :: plotVarLabels


  if (ubound(io_unklabelsGlobal,1) > MAX_MESH_VAR) then
     !This check is needed to prevent overflowing many statically size arrays.
     call Driver_abortFlash("Mesh replication support with type based I/O "//&
          "is currently limited.  The I/O code has statically sized arrays "//&
          "that will overflow.  You can make these arrays bigger and get "//&
          "past this abort by increasing MAX_MESH_VAR in io_flash.h")
  end if

  if (io_fileFormatVersion == 10) then
     call RuntimeParameters_get("meshCopyCount", meshCopyCount)
     if (meshCopyCount > 1) then
        call Driver_abortFlash("Mesh replication not supported for fileFmt=10")
     end if
  end if


  !NEW file format
  !*****************************************************************************
  !*****************************************************************************
  eachStruct: do i = lbound(allGridDataStruct,1), ubound(allGridDataStruct,1)

     gridDataStruct = allGridDataStruct(i)

     !We assume that each block belonging to a single process has the same size.
     !However, we allow different processes to have different block sizes.
     call io_getZeroBasedBlkSubarray(gridDataStruct, blockInnerSize, &
          blockOuterSize, blockInnerOffset)

     call io_getZeroBasedVarInfo(PLOTFILE, gridDataStruct, numGridVars, &
          numPlotVars, plotVarOffsets, plotVarLabels)

     call Grid_getNumVars(gridDataStruct, extentVarDim)

     !The fileFmt=10 types are invalid when there is mesh replication.
     call io_init_grid_mpi_types(myPE, gridDataStruct, blockOuterSize, &
          blockInnerSize, blockInnerOffset, extentVarDim, &
          plotVarOffsets, numPlotVars)

  end do eachStruct
  !Need the MPI types created in io_init_grid_mpi_types for restarting.
  !*****************************************************************************
  !*****************************************************************************

#ifdef FLASH_IO_HDF5
  call RuntimeParameters_get('packMeshPlotWriteHDF5', io_packMeshPlotWriteHDF5)
  call RuntimeParameters_get('packMeshChkWriteHDF5', io_packMeshChkWriteHDF5)
  call RuntimeParameters_get('packMeshChkReadHDF5', io_packMeshChkReadHDF5)
#else
  io_packMeshPlotWriteHDF5 = .false.
  io_packMeshChkWriteHDF5 = .false.
  io_packMeshChkReadHDF5 = .false.
#endif

#ifdef FLASH_IO_PNETCDF
  call RuntimeParameters_get("asyncMeshPlotWritePnet", io_asyncMeshPlotWritePnet)
  call RuntimeParameters_get("asyncMeshChkWritePnet", io_asyncMeshChkWritePnet)
  call RuntimeParameters_get("asyncMeshChkReadPnet", io_asyncMeshChkReadPnet)

  !If any of the above 3 pnetcdf runtime parameters are true then
  !we must initialize code specific to nonblocking pnetcdf feature.
  if (io_asyncMeshPlotWritePnet .or. io_asyncMeshChkWritePnet .or. &
       io_asyncMeshChkReadPnet) then
     call io_ncmpi_nonblocking_init(myPE)
  end if
#else
  io_asyncMeshPlotWritePnet = .false.
  io_asyncMeshChkWritePnet = .false.
  io_asyncMeshChkReadPnet = .false.
#endif

  call RuntimeParameters_get("useLegacyLabels", io_useLegacyLabels)

end subroutine io_typeInit
