!!****if* source/IO/IOTypes/io_xfer_mesh_data
!!
!! NAME
!!  io_xfer_mesh_data
!!
!! SYNOPSIS
!!
!!  io_xfer_mesh_data(integer(in) :: fileID,
!!                    integer(in) :: fileFmt,
!!                    integer(in) :: fileType,
!!                    integer(in) :: libType,
!!                    integer(in) :: xferType,
!!                    integer(in) :: localNumBlocks,
!!                    integer(in) :: globalBlockOffset)
!!
!!
!! DESCRIPTION
!!
!! This subroutine calls a C function to transfer mesh data
!! between memory and file.  It works by setting a pointer to 
!! UNK, FACEX, FACEY, FACEZ and SCRATCH arrays, and then passing the
!! address of the first grid element to the C transfer function.
!! It can transfer data from memory to file and from file to memory.
!!
!!
!! ARGUMENTS
!!
!! fileID: The HDF5 or pnetcdf file identifier (used directly by the libraries)
!! fileFmt: The FLASH file layout.  The standard layout that tools
!!          understand is 9 (one mesh variable per dataset), but there is also
!!          support for an experimental file layout 10 (all mesh variables in
!!          the same dataset)
!! fileType: The FLASH file type (checkpoint file or plot file)
!! libType: The library we are using (HDF5 or pnetcdf)
!! xferType: The direction of data transfer: memory to file or file to memory
!! localNumBlocks: The number of blocks on myPE being trasferred to/from file
!! globalBlockOffset: The read/write block offset in file.
!!
!!***


#include "constants.h"
#include "Flash.h"
#include "io_flash.h"

subroutine io_xfer_mesh_data(fileID, fileFmt, fileType, &
     libType, xferType, localNumBlocks, globalBlockOffset)

#ifdef USE_IO_C_INTERFACE
  use iso_c_binding, ONLY : c_loc, c_char
  use io_c_interface, ONLY : io_attribute_write, &
       io_ncmpi_define_mode_enddef, &
       io_ncmpi_define_mode_redef
  use io_c_type_interface, ONLY : io_ncmpi_nonblocking_complete_requests, &
       io_xfer_mesh_dataset
#endif
  use physicaldata, ONLY : unk, facevarx, facevary, facevarz
  use Grid_data, ONLY : scratch
  use io_typeInterface, ONLY : io_getZeroBasedBlkSubarray, &
       io_getZeroBasedVarInfo, io_do_xfer
  use io_typeData, ONLY : io_packMeshPlotWriteHDF5, io_packMeshChkWriteHDF5, &
       io_packMeshChkReadHDF5, io_asyncMeshPlotWritePnet, &
       io_asyncMeshChkWritePnet, io_asyncMeshChkReadPnet
  use Grid_interface, ONLY : Grid_getBlkCornerID, Grid_getNumVars
  use Driver_interface, ONLY : Driver_abortFlash
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Simulation_interface, ONLY : Simulation_mapStrToInt
  use IO_data, ONLY : io_globalMe

  implicit none
  integer, intent(IN) :: fileID, fileFmt, fileType, libType, &
       xferType, localNumBlocks, globalBlockOffset
#ifndef USE_IO_C_INTERFACE
#define c_loc(x) x
  integer, parameter :: c_char = KIND('A')
#endif

  real, pointer, dimension(:,:,:,:,:) :: gridData
  real, dimension(MAX_MESH_VAR) :: globalVarMinAll, globalVarMaxAll
  real, dimension(MAX_MESH_VAR), target :: globalVarMinSubset, &
       globalVarMaxSubset !target required for c_loc.

  integer, dimension(MAX_MESH_VAR) :: globalVarOffsets, localVarOffsets
  integer, dimension(MDIM) :: blockInnerSize, blockOuterSize, blockInnerOffset, &
       globalIndexLimits
  integer, dimension(IO_MAX_DIMS) :: array_of_sizes, array_of_subsizes, &
       array_of_subsizes_arg, array_of_starts, globalOffsetArray
  integer :: numMeshVar, numXferMeshVar, d, i, j, varIndex

  integer, parameter, dimension(5) :: gridDataStructs = &
       (/CENTER, SCRATCH, FACEX, FACEY, FACEZ/)

  character (kind=c_char,len=MAX_STRING_LENGTH), dimension(MAX_MESH_VAR) :: &
       meshVarLabels
  character (kind=c_char,len=MAX_STRING_LENGTH) :: dset_str, xfer_str

  character(kind=c_char,len=7), parameter :: mesh_str_tbl(5) = (/ &
       'unknown', &
       'facex  ', &
       'facey  ', &
       'facez  ', &
       'scratch' /)
  character(kind=c_char,len=7) :: min_att_str = "minimum"
  character(kind=c_char,len=7) :: max_att_str = "maximum"
  integer :: nonBlocking, numFileDatasets, prePackData, numDataDims
  integer :: cornerID(MDIM), stride(MDIM), extentVarDim
  logical :: doXfer
#ifdef DEBUG_IO
  logical, parameter :: debugIO = .true.
#else
  logical, parameter :: debugIO = .false.
#endif

  if (xferType == IO_WRITE_XFER) then
     xfer_str = "write mesh data"
  else if (xferType == IO_READ_XFER) then
     xfer_str = "read mesh data"
  else
     call Driver_abortFlash("Invalid transfer type")
  end if
  IO_TIMERS_START(xfer_str)
  nullify(gridData)

  !Do we use non blocking transfers (currently applicable to Pnetcdf only).
  nonBlocking = 0
  if (libType == IO_FILE_PNETCDF) then
     if (fileType == CHECKPOINTFILE) then 
        if ( (xferType == IO_READ_XFER .and. io_asyncMeshChkReadPnet) .or. &
             (xferType == IO_WRITE_XFER .and. io_asyncMeshChkWritePnet)) then
           nonBlocking = 1
        end if
     else if (fileType == PLOTFILE) then 
        if ( (xferType == IO_WRITE_XFER .and. io_asyncMeshPlotWritePnet)) then
           nonBlocking = 1
        end if
     end if
  end if


  !Do we pack mesh data (currently applicable to HDF5 only).
  prePackData = 0
  if (libType == IO_FILE_HDF5) then
     if (fileType == CHECKPOINTFILE) then 
        if ( (xferType == IO_READ_XFER .and. io_packMeshChkReadHDF5) .or. &
             (xferType == IO_WRITE_XFER .and. io_packMeshChkWriteHDF5)) then
           prePackData = 1
        end if
     else if (fileType == PLOTFILE) then 
        if ( (xferType == IO_WRITE_XFER .and. io_packMeshPlotWriteHDF5)) then
           prePackData = 1
        end if        
     end if
  end if

  if (debugIO) then
     if (io_globalMe == MASTER_PE) then
        write(6,'(a,i2,a,i2)') &
             " [io_xfer_mesh_data]: *** Using non-blocking writes:", &
             nonBlocking, ", pre-packing data:", prePackData
     end if
  end if

  eachDataStruct: do d = lbound(gridDataStructs,1), ubound(gridDataStructs,1)
     select case (gridDataStructs(d))
     case (CENTER)
        gridData => unk
     case (FACEX)
        gridData => facevarx
     case (FACEY)
        if (NDIM > 1) then
           gridData => facevary
        else
           nullify(gridData) !So we crash early if we mistakenly use it.
        end if
     case (FACEZ)
        if (NDIM > 2) then
           gridData => facevarz
        else
           nullify(gridData) !So we crash early if we mistakenly use it.
        end if
     case (SCRATCH)
        gridData => scratch
     case DEFAULT
        call Driver_abortFlash("Data structure not recognised")
     end select

     call Grid_getNumVars(gridDataStructs(d), extentVarDim)

     !Get global info about all mesh variables to be transferred to/from file.
     call io_getZeroBasedVarInfo(fileType, gridDataStructs(d), numMeshVar, &
          numXferMeshVar, globalVarOffsets, meshVarLabels)

     !Now that we support mesh replication (CENTER only) the global offsets
     !may not be the same as the local offsets.
     if (gridDataStructs(d) == CENTER) then
        do i = 1, numXferMeshVar
           call Simulation_mapStrToInt(trim(meshVarLabels(i)), &
                varIndex, MAPBLOCK_UNK)
           if (varIndex /= NONEXISTENT) then
              localVarOffsets(i) = varIndex - 1
           else
              localVarOffsets(i) = NONEXISTENT
           end if
        end do
     else
        localVarOffsets(1:numXferMeshVar) = globalVarOffsets(1:numXferMeshVar)
     end if

     !WARNING. numMeshVar and numXferMeshVar are obtained in a global view.
     !This means we cannot assume that they are correct in the local view
     !when using mesh replication.  Note that they are passed (directly and
     !indirectly) to io_xfer_mesh_dataset, but this does not cause a
     !problem because they are only used by the lower-level routine
     !io_h5_xfer_mesh_dataset for fileFmt=10.  We do not support mesh
     !replication for fileFmt=10 (and we should have already aborted in
     !io_typeInit if we tried) so I do not mind passing these
     !misleading values.  The code needs to be refactored with the
     !requirements of mesh replication in mind.

     CheckForFileXfer: if (numXferMeshVar > 0) then
        if (debugIO) then
           if (io_globalMe == MASTER_PE) then
              if (fileType == CHECKPOINTFILE) &
                   print *, "[io_xfer_mesh_data]: Checkpoint file transfer of ", &
                   trim(mesh_str_tbl(d)), " with file format ", fileFmt
              if (fileType == PLOTFILE) &
                   print *, "[io_xfer_mesh_data]: Plot file transfer of ", &
                   trim(mesh_str_tbl(d)), " with file format ", fileFmt
      
              do i = 1, numXferMeshVar
                 print *, "[io_xfer_mesh_data]: Will transfer variable ", &
                      trim(meshVarLabels(i)), " at memory offset", globalVarOffsets(i)
              end do
           end if
        end if

        !We assume that each block belonging to a single process has the same size.
        !However, we allow different processes to have different block sizes.
        call io_getZeroBasedBlkSubarray(gridDataStructs(d), blockInnerSize, &
             blockOuterSize, blockInnerOffset)
        
        !Total amount of local data as stored in memory:
        array_of_sizes(1) = localNumBlocks
        array_of_sizes(2) = blockOuterSize(3)
        array_of_sizes(3) = blockOuterSize(2)
        array_of_sizes(4) = blockOuterSize(1)
        array_of_sizes(5) = extentVarDim
        
        !Total amount of local data to transfer:
        array_of_subsizes(1) = localNumBlocks
        array_of_subsizes(2) = blockInnerSize(3)
        array_of_subsizes(3) = blockInnerSize(2)
        array_of_subsizes(4) = blockInnerSize(1)
        array_of_subsizes(5) = numXferMeshVar
        
        !Local offset in memory to reach data to transfer:
        array_of_starts(1) = 0
        array_of_starts(2) = blockInnerOffset(3)
        array_of_starts(3) = blockInnerOffset(2)
        array_of_starts(4) = blockInnerOffset(1)
        array_of_starts(5) = 0
   
        !We need to be aware that NOFBS grid applications store data in file
        !differently to FBS grid applications (i.e. PM and UG).  In NOFBS mode
        !the grid data is stored in one large block.
#ifdef IO_FLASH_NOFBS_UG
        !Find my file offset.  Always zero unless X,Y,Z dimension.
        globalOffsetArray(1) = 0
        call Grid_getBlkCornerID(1, cornerID, stride)
        globalOffsetArray(2) = (cornerID(KAXIS)-1) * K3D
        globalOffsetArray(3) = (cornerID(JAXIS)-1) * K2D
        globalOffsetArray(4) = cornerID(IAXIS)-1
        globalOffsetArray(5) = 0
#else
        !Find my file offset.  Always zero unless block dimension.
        globalOffsetArray(1) = globalBlockOffset
        globalOffsetArray(2:5) = 0
#endif


        !To help spot errors I fill the guardcells with -8.
        !unk(:,1 : NGUARD,:,:,:) = -8.0
        !unk(:,NXB+NGUARD+1 : NXB+(2*NGUARD),:,:,:) = -8.0
        !if(NDIM >= 2) then
        !   unk(:,:,1 : K2D*NGUARD,:,:) = -8.0
        !   unk(:,:,NYB+K2D*(NGUARD+1) : NYB+K2D*(2*NGUARD),:,:) = -8.0
        !end if
        !if(NDIM == 3) then
        !   unk(:,:,:,1 : K3D*NGUARD,:) = -8.0
        !   unk(:,:,:,NZB+K3D*(NGUARD+1) : NZB+K3D*(2*NGUARD),:) = -8.0
        !end if



        !Get the max and minimum values (NOTE: Modify this routine so we don't 
        !generate values for grid variables that are not output).
        if (xferType == IO_WRITE_XFER) then
           if (numXferMeshVar > 0) then
              call io_getVarExtrema(numMeshVar, globalVarMinAll, globalVarMaxAll, &
                   gridDataStructs(d))
              do i = 1, numXferMeshVar
                 globalVarMinSubset(i) = globalVarMinAll(globalVarOffsets(i) + 1)
                 globalVarMaxSubset(i) = globalVarMaxAll(globalVarOffsets(i) + 1)
              end do
           end if
        end if



        if (fileFmt <= 9) then
           !i.e. N read/writes of 1 variable.
           numFileDatasets = numXferMeshVar
           numDataDims = 4
        else if (fileFmt == 10) then
           !i.e. 1 read/write of N variables.
           numFileDatasets = 1
           numDataDims = 5
        else
           call Driver_abortFlash("File format not recognised")
        end if


        do i = 1, numFileDatasets

           !We must always pass a valid memory address from gridData (hence j=1)
           j = 1
           array_of_subsizes_arg = array_of_subsizes


           if (fileFmt <= 9) then
              call io_do_xfer(xferType, gridDataStructs(d), meshVarLabels(i), &
                   doXfer)

              dset_str = trim(meshVarLabels(i))
              if (doXfer) then
                 j = localVarOffsets(i) + 1 !gridData element (unit-based).
                 if (j < lbound(gridData,1) .or. j > ubound(gridData,1)) then
                    call Driver_abortFlash("Local index out of bounds")
                 end if
                 if (debugIO) then
                    if (io_globalMe == MASTER_PE) then
                       write(6,'(a,i4,a,a,a,i4)')" Processor ", io_globalMe, &
                            " will transfer ", trim(dset_str), " to/from ", j
                    end if
                 end if
              else
                 !Setting the subsizes to zero prevents my C code from
                 !transferring data.
                 array_of_subsizes_arg = 0
              end if

           else if (fileFmt == 10) then
              !Mesh replication does not work with fileFmt=10.  There is no
              !error checking here because we should have already aborted
              !in io_typeInit.  WARNING!
              dset_str = trim(mesh_str_tbl(d))
           end if


           call io_xfer_mesh_dataset(&
                io_globalMe, &
                fileID, &
                libType, &
                xferType, &
                fileType, &
                fileFmt, &
                gridDataStructs(d), &
                numDataDims, &
                numXferMeshVar, &
                nonBlocking, &
                prePackData, &
                globalOffsetArray, &
                array_of_sizes, &
                array_of_subsizes_arg, &
                array_of_starts, &
                localVarOffsets, &
                dset_str, &
                len_trim(dset_str), &
                c_loc(gridData(j,1,1,1,1)))


           if (xferType == IO_WRITE_XFER) then

#if defined(FLASH_IO_PNETCDF)
              call io_ncmpi_define_mode_redef(fileID)
#endif

              call io_attribute_write(io_globalMe, fileID, libType, &
                   IO_FLASH_DOUBLE, dset_str, len_trim(dset_str), &
                   min_att_str, len_trim(min_att_str), &
                   c_loc(globalVarMinSubset(i)))
              call io_attribute_write(io_globalMe, fileID, libType, &
                   IO_FLASH_DOUBLE, dset_str, len_trim(dset_str), &
                   max_att_str, len_trim(max_att_str), &
                   c_loc(globalVarMaxSubset(i)))

#if defined(FLASH_IO_PNETCDF)
              call io_ncmpi_define_mode_enddef(fileID)
#endif

           end if
        end do

        if (nonBlocking == 1) then

#if defined(FLASH_IO_PNETCDF)
           call io_ncmpi_nonblocking_complete_requests(fileID);
#endif

        end if

     endif CheckForFileXfer
  end do eachDataStruct

  IO_TIMERS_STOP(xfer_str)
end subroutine io_xfer_mesh_data
