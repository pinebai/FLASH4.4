!!****if* source/IO/IOMain/hdf5/parallel/PM_argonne/io_writeData
!!
!! NAME
!!
!!  io_writeData
!!
!!
!! SYNOPSIS
!!
!!  io_writeData(integer(in) :: fileID) 
!!          
!!
!!
!!
!! DESCRIPTION
!!
!!  This function writes the checkpoint data to an hdf5 file to store the 
!!  paramesh data.  IO is done in parallel -- no copying of the data to 
!!  a single processor
!!  to do the writing is performed.  HDF5 v. 1.4.0 or later is required
!!
!!  HDF5 uses MPI-IO (via ROMIO) to support parallel IO.  Each processor
!!  must open the file, create the datasets, and dataspaces for each HDF
!!  record.
!!
!!  A single record for each of the PARAMESH data structures is created.  A
!!  processor only writes to a subset of this record.  Each record has a
!!  dimension with length = globalNumBlocks.  The offset of a processor into this
!!  dimension is computed by looking at the total number of blocks that are
!!  below the current processor.
!!
!!  In this version of io_writeData, each variable is given its own
!!  record -- this makes it easier to change the variable list in the
!!  future without disturbing the format of the file.
!!
!!
!! ARGUMENTS
!! 
!!  fileID - integer file identifier for hdf5 file
!!
!!
!! NOTES
!!  variables that start with "io_" belong to the data module IO_data.
!!  the "io_" is meant to indicate that these variables have IO unit 
!!  scope.   For performance purposes IO particularly, io_writeData
!!  and read data are allowed to access unk directly.  Additionally,
!!  io_writeData needs access to some variables that are specific to 
!!  Grid.  These are stored in Grid_ioData and variables associated
!!  with Grid_ioData start with "gio_"
!!
!!***

#include "constants.h"
#include "Flash.h"
#include "io_flash.h"

subroutine io_writeData(fileID) 

#ifdef USE_IO_C_INTERFACE
  use iso_c_binding, ONLY : c_loc, c_char
#endif

  use IO_data, ONLY : io_globalMe, io_globalNumProcs,  io_realParmNames, io_realParmValues, io_numRealParms, &
       io_intParmNames, io_intParmValues, io_numIntParms, &
       io_logParmNames, io_logParmValues, io_numLogParms, &
       io_strParmNames, io_strParmValues, io_numStrParms, &
       io_realScalarNames, io_realScalarValues, io_numRealScalars, &
       io_intScalarNames, io_intScalarValues, io_numIntScalars, &
       io_logScalarNames, io_logScalarValues, io_numLogScalars, &
       io_strScalarNames, io_strScalarValues, io_numStrScalars, &
       io_logToIntScalarValues, io_logToIntParmValues, io_unklabels, &
       io_geometry, io_setupCall, io_buildDir, io_flashRelease, &
       io_fileCreationTime, io_buildDate, io_buildMachine, io_cflags, io_fflags, &
       io_setupTimeStamp, io_buildTimeStamp, io_doublePrecision, &
       io_ilo, io_ihi, io_jlo, io_jhi, io_klo, io_khi, &
       io_plotVarStr, io_nPlotVars, io_outputSplitNum, &
       io_chkGuardCellsInput, io_chkGuardCellsOutput, &
       io_plotGridVarStr, io_faceXVarLabels, io_faceYVarLabels, &
       io_faceZVarLabels, io_comm, io_splitNumBlks, io_useCollectiveHDF5, &
       io_plotfileGridQuantityDP, io_plotfileMetadataDP, io_fileFormatVersion, &
       tree_data_t, io_meshMe, io_acrossMe

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_abortFlash
  use Logfile_interface, ONLY : Logfile_stamp
  use Grid_interface, ONLY : Grid_getLocalNumBlks, Grid_getBlkIndexLimits, &
       Grid_getGlobalIndexLimits, Grid_getBlkBoundBox, &
       Grid_getBlkCenterCoords, Grid_getBlkPhysicalSize
#ifdef FLASH_GRID_PARAMESH
  use tree, ONLY : nodetype, bnd_box, coord, bsize, lrefine, nfaces, nchild
#ifdef FLASH_GRID_PARAMESH3OR4
  use tree, ONLY : MFLAGS, which_child, bflags
  use Grid_data, ONLY : gr_gsurr_blks
#endif
#endif
  use IO_interface, ONLY : IO_setScalar
  use Grid_data, ONLY : gr_globalNumBlocks, gr_nToLeft, gr_gid, &
       gr_globalOffset, scratch
  use io_typeInterface, ONLY : io_getZeroBasedBlkSubarray, &
       io_getZeroBasedVarInfo, io_xfer_mesh_data, io_xfer_tree_data, &
       io_create_grid_header
  use Timers_interface, ONLY : Timers_start, Timers_stop
  implicit none


#include "Flash_mpi.h"

  integer, intent(in) :: fileID
  type (tree_data_t) :: tree_data

#ifndef USE_IO_C_INTERFACE
#define c_loc(x) x
  integer, parameter :: c_char = KIND('A')
#endif

  character (len=MAX_STRING_LENGTH) :: buff
  character(len=16) :: numToStr  

  !Adjust our offset so that we can eliminate redundant data in a group
  integer :: localNumBlocks, localOffset, splitOffset, ierr

  real, dimension(2,MDIM) :: boundBox
  real, dimension(MDIM) :: blockCenterCoords, blockSize

  !presentDims needs to depend on FILE_FORMAT_VERSION.
  integer, parameter :: xferType = IO_WRITE_XFER, presentDims = MDIM, &
       libType = IO_FILE_HDF5
  integer :: dataFP, attributeFP, fileType, fileFmt, numFileBlks

#ifndef FLASH_GRID_PARAMESH
  !Declare arrays that exist only in Paramesh simulations.
  integer, target, dimension(1) :: nodetype = (/LEAF/), lrefine = (/1/)
  real, target, dimension(2,MDIM,1) :: bnd_box
  real, target, dimension(MDIM,1) :: coord, bsize
#endif
#ifdef IO_FLASH_NOFBS_UG
  integer, dimension(MDIM) :: globalIndexLimits

  !get the global index limits of the domain and set the
  !values of nxb, nyb and nzb in the scalar list.
  call Grid_getGlobalIndexLimits(globalIndexLimits)
  call IO_setScalar("nxb", globalIndexLimits(IAXIS))
  call IO_setScalar("nyb", globalIndexLimits(JAXIS))
  call IO_setScalar("nzb", globalIndexLimits(KAXIS))
  call IO_setScalar("globalNumBlocks", 1)
#endif

  IO_TIMERS_START("io_writeData")
  fileFmt = io_fileFormatVersion

  !Floating point type for mesh datasets + mesh dataset attributes
  if (io_doublePrecision .or. io_plotfileGridQuantityDP) then
     dataFP = IO_FLASH_DOUBLE
     attributeFP = IO_FLASH_DOUBLE
  else
     dataFP = IO_FLASH_FLOAT
     attributeFP = IO_FLASH_FLOAT
  end if



  !! call the generic function prepareLists to allocate and 
  !! fill the runtime parameter lists and the scalar lists
  IO_TIMERS_START("write lists/sim info")
  call io_prepareListsWrite()


  call Grid_getLocalNumBlks(localNumBlocks)


  !Find our local offset for split-file IO
  if (io_outputSplitNum > 1) then

     call MPI_ALLREDUCE(gr_globalOffset, splitOffset, 1, FLASH_INTEGER, &
          MPI_MIN, io_comm, ierr)
     localOffset = gr_globalOffset - splitOffset
     !find number of blocks for a file
     call MPI_ALLREDUCE(localNumBlocks, io_splitNumBlks, 1, FLASH_INTEGER,&
          MPI_SUM, io_comm, ierr)
  else
     localOffset = gr_globalOffset
     io_splitNumBlks = gr_globalNumBlocks
  end if

#if defined(IO_FLASH_NOFBS_UG)
  numFileBlks = 1
#else
  numFileBlks = io_splitNumBlks
#endif

  !If we are writing a checkpoint file (io_doublePrecision == .true.)
  !!write all of the unk vars.  
  !If writing a plotfile then only certain variables are written
  if(io_doublePrecision) then

     !! write the header info
     call io_h5write_header(io_globalMe, NUNK_VARS, fileID, io_geometry, &
          io_unklabels, io_setupCall, io_fileCreationTime, io_flashRelease, &
          io_buildDate, io_buildDir, io_buildMachine, io_cflags, io_fflags, &
          io_setupTimeStamp, io_buildTimeStamp, fileFmt, io_outputSplitNum)

  else

     call io_h5write_header(io_globalMe, io_nPlotVars, fileID, io_geometry, &
          io_plotVarStr, io_setupCall, io_fileCreationTime, io_flashRelease, &
          io_buildDate, io_buildDir, io_buildMachine, io_cflags, io_fflags, &
          io_setupTimeStamp, io_buildTimeStamp, fileFmt, io_outputSplitNum)

  end if


  !! write the runtime parameters
  call io_h5write_runtime_parameters(io_globalMe, &
       fileID, &
       io_numRealParms, &
       io_realParmNames, &
       io_realParmValues, &
       io_numIntParms, &
       io_intParmNames, &
       io_intParmValues, &
       io_numStrParms, &
       io_strParmNames, &
       io_strParmValues, &
       io_numLogParms, &
       io_logParmNames, &
       io_logToIntParmValues, &
       io_outputSplitNum)



  !! write the scalars
  call io_h5write_scalars(io_globalMe, &
       fileID, &
       io_numRealScalars, &
       io_realScalarNames, &
       io_realScalarValues, &
       io_numIntScalars, &
       io_intScalarNames, &
       io_intScalarValues, &
       io_numStrScalars, &
       io_strScalarNames, &
       io_strScalarValues, &
       io_numLogScalars, &
       io_logScalarNames, &
       io_logToIntScalarValues, &
       io_outputSplitNum)


  call io_finalizeListsWrite()
  IO_TIMERS_STOP("write lists/sim info")


  if(io_doublePrecision) then  
     fileType = CHECKPOINTFILE
  else
     fileType = PLOTFILE
  end if


  !! Create all the datasets up front.  This allows us to use the same
  !! data transfer interface as pnetcdf.
  IO_TIMERS_START("create tree datasets")
  call io_createDatasets(fileID, numFileBlks, presentDims)
  IO_TIMERS_STOP("create tree datasets")


  IO_TIMERS_START("create mesh datasets")
  call io_create_grid_header(io_globalMe, fileID, fileFmt, fileType, &
       libType, dataFP, attributeFP)
  IO_TIMERS_STOP("create mesh datasets")


#ifdef IO_FLASH_NOFBS_UG
  !Metadata is only written by the master processor for NOFBS UG.
  !See r13578 NoFbs io_writeData for mesh replication logic.
  if (io_meshMe == MASTER_PE .and. io_acrossMe == 0) then
     localNumBlocks = 1
  else
     localNumBlocks = 0
  end if

  !if faking one block, everything is a boundary, no neighbors
  gr_gid = -21

  !use entire physical domain to 'fake' single block
  call RuntimeParameters_get("xmin", bnd_box(LOW,IAXIS,1))
  call RuntimeParameters_get("xmax", bnd_box(HIGH,IAXIS,1))
  call RuntimeParameters_get("ymin", bnd_box(LOW,JAXIS,1))
  call RuntimeParameters_get("ymax", bnd_box(HIGH,JAXIS,1))
  call RuntimeParameters_get("zmin", bnd_box(LOW,KAXIS,1))
  call RuntimeParameters_get("zmax", bnd_box(HIGH,KAXIS,1))
  coord(:,1) = bnd_box(HIGH,:,1) * 0.5
  bsize(:,1) = bnd_box(HIGH,:,1) - bnd_box(LOW,:,1)
#endif

#if defined(IO_FLASH_UG)
  call Grid_getBlkBoundBox(1, boundBox)
  bnd_box(:,:,1) = boundBox(:,:)

  call Grid_getBlkCenterCoords(1, blockCenterCoords)
  coord(:,1) = blockCenterCoords(:)

  call Grid_getBlkPhysicalSize(1, blockSize)
  bsize(:,1) = blockSize(:)
#endif

  tree_data % bnd_box => bnd_box
  tree_data % coord => coord
  tree_data % bsize => bsize
  tree_data % gid => gr_gid
  tree_data % nodetype => nodetype
  tree_data % lrefine => lrefine
#ifdef FLASH_GRID_PARAMESH3OR4
  tree_data % bflags => bflags
  tree_data % which_child => which_child
  tree_data % gsurr_blks => gr_gsurr_blks
#else
  nullify(tree_data % bflags)
  nullify(tree_data % which_child)
  nullify(tree_data % gsurr_blks)
#endif
  allocate(tree_data % procnumber(max(1,localNumBlocks)))
  tree_data % procnumber(:) = io_meshMe

  call io_xfer_tree_data(tree_data, fileID, IO_FILE_HDF5, IO_WRITE_XFER, &
       localNumBlocks, localOffset, presentDims)

  deallocate(tree_data % procnumber)
  nullify(tree_data % procnumber)


#ifdef IO_FLASH_NOFBS_UG
  localNumBlocks = 1
#endif

  call io_xfer_mesh_data(fileID, fileFmt, fileType, IO_FILE_HDF5, &
       IO_WRITE_XFER, localNumBlocks, localOffset)


  if (io_globalMe .EQ. MASTER_PE) then
     write (numToStr, "(I7)") gr_globalNumBlocks
     buff = "wrote " // numToStr // " blocks"
     call Logfile_stamp( buff, '[io_writeData]')
  end if

  IO_TIMERS_STOP("io_writeData")
end subroutine io_writeData
