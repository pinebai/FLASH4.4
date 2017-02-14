!!****if* source/IO/IOMain/pnetcdf/io_writeData
!!
!! NAME
!!
!!  io_writeData
!!
!!
!! SYNOPSIS
!!
!!  io_writeData() 
!!                
!!               integer(in) :: fileID)
!!
!!
!! DESCRIPTION
!!
!!
!!  This version of io_writeData uses Parallel Netcdf to store the 
!!  data.  The 
!!  IO is done in parallel -- no copying of the data to a single processor
!!  to do the writing is performed.  Pnetcdf v. 1.0.0 or later is required
!!
!!  Pnetcdf uses MPI-IO (via ROMIO) to support parallel IO.  Each processor
!!  must open the file, create the datasets, and dataspaces for each netcdf 
!!  record.  
!!
!!  A single record for each data structure is created.  A
!!  processor only writes to a subset of this record.  Each record has a 
!!  dimension with length = tot_blocks.  The offset of a processor into this 
!!  dimension is computed by looking at the total number of blocks that are
!!  below the current processor.
!!
!!  In this version of the checkpoint, each variable is given its own 
!!  record -- this makes it easier to change the variable list in the
!!  future without disturbing the format of the file.  
!!
!!
!! ARGUMENTS
!! 
!!  fileID - integer file identifier for pnetcdf file
!!
!!***

#include "constants.h"
#include "Flash.h"
#include "io_flash.h"

subroutine io_writeData(fileID)

#ifdef USE_IO_C_INTERFACE
  use iso_c_binding, ONLY : c_loc, c_char
  use io_c_interface, ONLY : io_ncmpi_read_file_format, &
       io_ncmpi_define_mode_enddef
#endif

  use IO_data, ONLY : io_globalMe, io_realParmNames, io_realParmValues, io_numRealParms, &
       io_intParmNames, io_intParmValues, io_numIntParms, &
       io_logParmNames, io_logParmValues, io_numLogParms, &
       io_strParmNames, io_strParmValues, io_numStrParms, &
       io_realScalarNames, io_realScalarValues, io_numRealScalars, &
       io_intScalarNames, io_intScalarValues, io_numIntScalars, &
       io_logScalarNames, io_logScalarValues, io_numLogScalars, &
       io_strScalarNames, io_strScalarValues, io_numStrScalars, &
       io_logToIntScalarValues, io_logToIntParmValues, &
       io_ilo, io_ihi, io_jlo, io_jhi, io_klo, io_khi, io_unklabels, &
       io_geometry, io_setupCall, io_buildDir, io_flashRelease, &
       io_fileCreationTime, io_buildDate, io_buildMachine, io_cflags, io_fflags, &
       io_setupTimeStamp, io_buildTimeStamp, io_doublePrecision, io_plotVarStr, &
       io_nPlotVars, io_maxPlotVars, io_plotGridVarStr, io_faceXVarLabels, io_faceYVarLabels, io_faceZVarLabels,&
       io_plotfileMetadataDP, io_plotfileGridQuantityDP, io_plotVar, &
       io_fileFormatVersion, io_globalNumProcs, &
       tree_data_t, io_meshMe, io_acrossMe
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_abortFlash
  use Particles_interface, ONLY : Particles_getGlobalNum
  use Grid_interface, ONLY : Grid_getLocalNumBlks, Grid_getBlkIndexLimits, &
       Grid_getGlobalIndexLimits, Grid_getBlkBoundBox, &
       Grid_getBlkCenterCoords, Grid_getBlkPhysicalSize
#ifdef FLASH_GRID_PARAMESH
  use tree, ONLY : nodetype, bnd_box, coord, bsize, lrefine, nfaces, nchild
#ifdef FLASH_GRID_PARAMESH3OR4
  use tree, ONLY : which_child, bflags
  use Grid_data, ONLY : gr_gsurr_blks
#endif
#endif
  use IO_interface, ONLY : IO_setScalar
  use Grid_data, ONLY : gr_globalNumBlocks, gr_globalOffset, gr_gid, scratch

  use physicaldata, ONLY : unk, facevarx, facevary, facevarz
  use io_typeInterface, ONLY : io_getZeroBasedBlkSubarray, &
       io_getZeroBasedVarInfo, io_xfer_mesh_data, io_xfer_tree_data, &
       io_create_grid_header
  implicit none

#include "Flash_mpi.h"

  type(tree_data_t) :: tree_data
  integer, intent(in) :: fileID
#ifndef USE_IO_C_INTERFACE
#define c_loc(x) x
  integer, parameter :: c_char = KIND('A')
#endif

  integer :: localNumBlocks, ngid, blockID, nOut
  integer :: i,j,lb, doublePrecisionInt, globalNumParticles
  integer :: nzonesBlock(MDIM)

  logical :: includeVar
  character (len=4), allocatable :: labels(:)
  integer :: metadataDPInt, gridQuantityDPInt

  integer, parameter :: libType = IO_FILE_PNETCDF, &
       xferType = IO_WRITE_XFER, presentDims = MDIM
  integer :: dataFP, attributeFP
  integer :: fileType, d
  integer :: fileFmt, totalBlocks
  integer :: numProcs

#ifndef FLASH_GRID_PARAMESH
  !Declare arrays that exist only in Paramesh simulations.
  integer, parameter :: nfaces = 2*NDIM, nchild = 2**NDIM
  integer, target, dimension(1) :: nodetype = (/LEAF/), lrefine = (/1/)
  real, target, dimension(2,MDIM,1) :: bnd_box
  real, target, dimension(MDIM,1) :: coord, bsize
  real, dimension(2,MDIM) :: boundBox
  real, dimension(MDIM) :: blockCenterCoords, blockSize
#endif

  numProcs = io_globalNumProcs
  fileFmt = io_fileFormatVersion


#ifdef IO_FLASH_NOFBS_UG
  !get the global index limits of the domain and set the
  !values of nxb, nyb and nzb in the scalar list.
  !nzonesBlock is used by io_ncmpi_write_header to set the
  !size of the pnetcdf variable on disk.  In NOFBS mode the
  !mesh array (stored as a pnetcdf variable) is faked up to be
  !1 block irrepective of the number of processors.
  call Grid_getGlobalIndexLimits(nzonesBlock)
  call IO_setScalar("nxb", nzonesBlock(IAXIS))
  call IO_setScalar("nyb", nzonesBlock(JAXIS))
  call IO_setScalar("nzb", nzonesBlock(KAXIS))
  call IO_setScalar("globalNumBlocks", 1)
  totalBlocks = 1
#else
  nzonesBlock(1) = NXB
  nzonesBlock(2) = NYB
  nzonesBlock(3) = NZB
  totalBlocks = gr_globalNumBlocks
#endif


  call io_prepareListsWrite()

  call Grid_getLocalNumBlks(localNumBlocks)

  if (NDIM == 1) then
     ngid = 5
  else if (NDIM == 2) then
     ngid = 9
  else if (NDIM == 3) then
     ngid = 15
  end if

  !are we writing a plotfile or checkpoint file?  make int equiv for c routine
  if(io_doublePrecision) then
     doublePrecisionInt = 1
     metadataDPInt = 1
     gridQuantityDPInt = 1
     nOut = NUNK_VARS
     allocate(labels(NUNK_VARS))
     labels = io_unklabels
     fileType = CHECKPOINTFILE
  else
     doublePrecisionInt = 0
     if (io_plotfileGridQuantityDP) then
        gridQuantityDPInt = 1
     else
        gridQuantityDPInt = 0
     endif
     if(io_plotfileMetadataDP) then
        metadataDPInt = 1
     else 
        metadataDPInt = 0
     endif
     nOut = io_nPlotVars
     allocate(labels(io_maxPlotVars))
     labels = io_plotVarStr  !Stores "none" in each element if no plot vars.
     fileType = PLOTFILE
  end if


  !Particles are not yet implemented for Paramesh
  call Particles_getGlobalNum(globalNumParticles)



  call io_ncmpi_write_header(io_globalMe, fileFmt, nOut, nzonesBlock, &
       labels, fileID, totalBlocks, ngid, &
       io_numRealParms, io_realParmNames, io_realParmValues, &
       io_numIntParms, io_intParmNames, io_intParmValues, &
       io_numStrParms, io_strParmNames, io_strParmValues, &
       io_numLogParms, io_logParmNames, io_logToIntParmValues, &
       io_numRealScalars, io_realScalarNames, io_realScalarValues, &
       io_numIntScalars, io_intScalarNames, io_intScalarValues, &
       io_numStrScalars, io_strScalarNames, io_strScalarValues, &
       io_numLogScalars, io_logScalarNames, io_logToIntScalarValues, &
       io_setupCall, io_fileCreationTime, io_flashRelease, io_buildDate, &
       io_buildDir, io_buildMachine, io_cflags, io_fflags, &
       io_setupTimeStamp,io_buildTimeStamp, &
       doublePrecisionInt, metadataDPInt, gridQuantityDPInt, globalNumParticles)


#ifdef FLASH_IO_EXPERIMENTAL
   if (1 == gridQuantityDPInt) then
     dataFP = IO_FLASH_DOUBLE
     attributeFP = IO_FLASH_DOUBLE
  else
     dataFP = IO_FLASH_FLOAT
     attributeFP = IO_FLASH_FLOAT
  end if

  call io_create_grid_header(io_globalMe, fileID, fileFmt, fileType, libType, &
       dataFP, attributeFP)

  !I now end define mode here.  Definitions are added in both 
  !io_ncmpi_write_header and io_ncmpi_write_grid_header.
  call io_ncmpi_define_mode_enddef(fileID)
#endif


  call io_finalizeListsWrite()
  !pnetcdf stores variables with an integer id.
  !the runtime parameters id = 0 and the scalars id = 1
  !The other variables the developer can handle since sometimes
  !additional variables may want to be checkpointed
  !The important thing to remember is that the variable id 
  !must be the same when reading in the data. It would also
  !be nice to keep the the ids the same in the UG and Paramesh cases


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

  call io_xfer_tree_data(tree_data, fileID, IO_FILE_PNETCDF, xferType, &
       localNumBlocks, gr_globalOffset, presentDims)

  deallocate(tree_data % procnumber)
  nullify(tree_data % procnumber)


#ifdef IO_FLASH_NOFBS_UG
  localNumBlocks = 1
#endif

  call io_xfer_mesh_data(fileID, fileFmt, fileType, &
       libType, xferType, localNumBlocks, gr_globalOffset)

  deallocate(labels)
end subroutine io_writeData
