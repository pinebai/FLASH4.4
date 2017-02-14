!!****if* source/IO/IOMain/hdf5/parallel/PM/io_writeData
!!
!! NAME
!!
!!  io_writeData
!!
!!
!! SYNOPSIS
!!
!!  call io_writeData(integer(in) :: fileID) 
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

!!REORDER(5):unk,facevar[xyz],scratch

#include "constants.h"
#include "Flash.h"
#include "io_flash.h"

subroutine io_writeData( fileID) 

#ifdef USE_IO_C_INTERFACE
  use iso_c_binding, ONLY : c_loc, c_char
#endif
  
  use IO_data, ONLY : io_globalMe, io_globalNumProcs,io_globalComm,&
         io_realParmNames, io_realParmValues, io_numRealParms, &
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
       io_chkGuardCellsOutput, &
       io_plotGridVarStr, io_faceXVarLabels, io_faceYVarLabels, &
       io_faceZVarLabels, io_comm, io_splitNumBlks, io_plotfileMetadataDP, &
       io_plotfileGridQuantityDP, io_fileFormatVersion, &
       io_meshMe, io_acrossMe, io_acrossNumProcs, io_unklabelsGlobal, io_unkNonRep, &
       tree_data_t

  use Simulation_interface, ONLY : Simulation_mapStrToInt
  use Driver_interface, ONLY : Driver_abortFlash
  use Logfile_interface, ONLY : Logfile_stamp
  use Grid_interface, ONLY : Grid_getLocalNumBlks, Grid_getBlkIndexLimits
     
  
  use Grid_data, ONLY : gr_globalNumBlocks, gr_nToLeft, &
       gr_globalOffset, scratch, gr_gid

#ifdef FLASH_GRID_PARAMESH
  use tree, ONLY : nodetype, lrefine, bnd_box, coord, bsize
#ifdef FLASH_GRID_PARAMESH3OR4
  use tree, ONLY : which_child, bflags
  use Grid_data, ONLY : gr_gsurr_blks
#endif
#endif
  
  use physicaldata, ONLY : unk, facevarx, facevary, facevarz
  use io_typeInterface, ONLY : io_xfer_tree_data

  implicit none
#ifndef USE_IO_C_INTERFACE
#define c_loc(x) x
  integer, parameter :: c_char = KIND('A')
#endif

#include "Flash_mpi.h"



  integer, intent(in) ::  fileID

  type (tree_data_t) :: tree_data
  character (len=MAX_STRING_LENGTH) :: buff
  character(len=16) :: numToStr  

  integer :: localNumBlocks, blockID
  integer :: i, lb, j, u, AllocateStatus
  integer :: dowrite
  integer, allocatable :: procnumber(:) 

  ! allocate storage to hold a single variable information
  ! this should only be a small memory overhead
  integer, parameter :: single = SELECTED_REAL_KIND(p=6)

!  real,save :: unkBufGC(1,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC,MAXBLOCKS)
  real,allocatable :: unkBufGC(:,:,:,:,:)

  ! allocate storage to hold the coordinate information and bounding box
  ! information
  real (kind=single) :: tmpSingle(MDIM,MAXBLOCKS)
  real (kind=single) :: bndSingle(2,MDIM,MAXBLOCKS)
  real (kind=single) :: spMax, spMin

  logical :: isPlotVar
  logical, allocatable :: isPlotVars(:)
  real, allocatable :: globalVarMin(:), globalVarMax(:)
  real, allocatable :: globalFaceXMin(:), globalFaceXMax(:)
  real, allocatable :: globalFaceYMin(:), globalFaceYMax(:)
  real, allocatable :: globalFaceZMin(:), globalFaceZMax(:)
  
  real, allocatable :: unkBuf(:,:,:,:,:)
  real, allocatable :: faceXBuf(:,:,:,:,:)
  real, allocatable :: faceYBuf(:,:,:,:,:)
  real, allocatable :: faceZBuf(:,:,:,:,:)
  
  real (kind=single),allocatable :: unkt(:,:,:,:,:)

  logical :: writeGuardCells = .true.

  integer :: blkLimits(2,MDIM), blkLimitsGC(2, MDIM)
  integer :: blkLimitsFX(2,MDIM), blkLimitsGcFX(2, MDIM)
  integer :: blkLimitsFY(2,MDIM), blkLimitsGcFY(2, MDIM)
  integer :: blkLimitsFZ(2,MDIM), blkLimitsGcFZ(2, MDIM)

  !Adjust our offset so that we can eliminate redundant data in a group
  integer :: localOffset, localRank, splitOffset
  integer :: ierr

  !presentDims needs to depend on FILE_FORMAT_VERSION.
  integer, parameter :: xferType = IO_WRITE_XFER, presentDims = MDIM, &
       libType = IO_FILE_HDF5
  integer :: numFileBlks


  !! call the generic function prepareLists to allocate and 
  !! fill the runtime parameter lists and the scalar lists
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

  numFileBlks = io_splitNumBlks
  
  
  !If we are writing a checkpoint file (io_doublePrecision == .true.)
  !!write all of the unk vars.  
  !If writing a plotfile then only certain variables are written
  if(io_doublePrecision) then

     !! write the header info
     call io_h5write_header(io_meshMe, ubound(io_unklabelsGlobal,1), fileID, io_geometry, &
          io_unklabelsGlobal, io_setupCall, io_fileCreationTime, io_flashRelease, &
          io_buildDate, io_buildDir, io_buildMachine, io_cflags, io_fflags, &
          io_setupTimeStamp, io_buildTimeStamp, io_fileFormatVersion, &
          io_outputSplitNum)
     
  else
     call io_h5write_header(io_meshMe, io_nPlotVars, fileID, io_geometry, &
          io_plotVarStr, io_setupCall, io_fileCreationTime, io_flashRelease, &
          io_buildDate, io_buildDir, io_buildMachine, io_cflags, io_fflags, &
          io_setupTimeStamp, io_buildTimeStamp, io_fileFormatVersion, &
          io_outputSplitNum)

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


  call io_createDatasets(fileID, numFileBlks, presentDims)



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




  allocate(globalVarMin(ubound(io_unklabelsGlobal,1)))
  allocate(globalVarMax(ubound(io_unklabelsGlobal,1)))

  !get the max and minimum variables
  call io_getVarExtrema(ubound(io_unklabelsGlobal,1), globalVarMin, globalVarMax, CENTER)

  


  !store the unknowns
  if (io_doublePrecision .or. io_plotfileGridQuantityDP) then
     if (.not. io_chkGuardCellsOutput) &
          allocate(unkBuf(1, NXB, NYB, NZB, MAXBLOCKS))
  else
     allocate(unkt(1, NXB, NYB, NZB, MAXBLOCKS))
  end if

  if (io_doublePrecision .and. io_chkGuardCellsOutput) then ! blkLimits is not used otherwise
     blockID = 1
     call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)
#ifndef FL_NON_PERMANENT_GUARDCELLS
     blkLimits = blkLimitsGC
#endif
  end if

  !do i = UNK_VARS_BEGIN,UNK_VARS_END
  do u=1, ubound(io_unklabelsGlobal,1)
     call Simulation_mapStrToInt(io_unklabelsGlobal(u),i,MAPBLOCK_UNK)
     ! only write mesh replicated data from mesh 0
     dowrite = 0
     if(i /= NONEXISTENT) then
        if(localNumBlocks > 0 .and. (io_acrossMe .eq. 0 .or. io_unkNonRep(i) > 0)) dowrite = 1
     else
        i = 1 ! gotta give it something
     end if
     
     if(io_doublePrecision) then
        if(io_chkGuardCellsOutput) then
           allocate(unkBufGC(1,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC,MAXBLOCKS))
           unkBufGC(1,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
                blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS), &
                blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS), 1:MAXBLOCKS) = &
                unk(i,:,:,:,1:MAXBLOCKS) 

           !DEV: The above works for permanent guardcells, though it leaves the 
           !guardcell region filled with arbitrary data in non-permanent
           !guardcell mode.
           call io_h5write_unknowns(io_globalMe, &
                fileID, & 
                GRID_IHI_GC, & 
                GRID_JHI_GC, & 
                GRID_KHI_GC, & 
                unkBufGC, & 
                globalVarMin(u), &
                globalVarMax(u), &
                io_unklabelsGlobal(u), &
                localNumBlocks, &
                io_splitNumBlks,  & 
                localOffset, &
                dowrite)
           deallocate(unkBufGC)
        else
           unkBuf(1,1:NXB,1:NYB,1:NZB,1:MAXBLOCKS) = &
                unk(i,io_ilo:io_ihi, io_jlo:io_jhi, io_klo:io_khi,1:MAXBLOCKS) 
           call io_h5write_unknowns(io_globalMe, &
                fileID, & 
                NXB, & 
                NYB, & 
                NZB, & 
                unkBuf, & 
                globalVarMin(u), &
                globalVarMax(u), &
                io_unklabelsGlobal(u), &
                localNumBlocks, &
                io_splitNumBlks,  & 
                localOffset, &
                dowrite)
        end if
     else
        isPlotVar = any(io_unklabelsGlobal(u) == io_plotVarStr(1:io_nPlotVars))
        if (isPlotVar) then
           if (io_plotfileGridQuantityDP) then
              unkBuf(1,1:NXB,1:NYB,1:NZB,1:MAXBLOCKS) = &
                   unk(i,io_ilo:io_ihi, io_jlo:io_jhi, io_klo:io_khi,1:MAXBLOCKS) 
              call io_h5write_unknowns(io_globalMe, &
                   fileID, & 
                   NXB, & 
                   NYB, & 
                   NZB, & 
                   unkBuf, & 
                   globalVarMin(u), &
                   globalVarMax(u), &
                   io_unklabelsGlobal(u), &
                   localNumBlocks, &
                   io_splitNumBlks,  & 
                   localOffset, &
                   dowrite)
           else
              unkt(1,:,:,:,1:localNumBlocks) = real(unk(i, & 
                   io_ilo:io_ihi, & 
                   io_jlo:io_jhi, & 
                   io_klo:io_khi, &
                   1:localNumBlocks), & 
                   kind = single)
              spMax = real(globalVarMax(u), kind = single)
              spMin = real(globalVarMin(u), kind = single)
              call io_h5write_unknowns_sp(io_globalMe, &
                   fileID, & 
                   NXB,   & 
                   NYB,   & 
                   NZB, & 
                   spMin, &
                   spMax, &
                   unkt,          & 
                   io_unklabelsGlobal(u),  & 
                   localNumBlocks,  & 
                   io_splitNumBlks,  & 
                   localOffset, &
                   dowrite)
           end if
        end if
     endif
     
     call MPI_BARRIER(io_globalComm, ierr)
  end do

  if (io_doublePrecision .or. io_plotfileGridQuantityDP) then
     if (.not. io_chkGuardCellsOutput) deallocate(unkBuf)
  else
     deallocate(unkt)
  end if
  deallocate(globalVarMin)
  deallocate(globalVarMax)

  if (io_doublePrecision .or. io_plotfileGridQuantityDP) then
     allocate(unkBuf(NXB, NYB, NZB, 1, 1:localNumBlocks))
  else
     allocate(unkt(NXB, NYB, NZB, 1, 1:localNumBlocks))
  end if
  allocate(globalVarMin(NSCRATCH_GRID_VARS))
  allocate(globalVarMax(NSCRATCH_GRID_VARS))
  allocate(isPlotVars(NSCRATCH_GRID_VARS))

  do i = SCRATCH_GRID_VARS_BEGIN,SCRATCH_GRID_VARS_END
     call io_isPlotVar(i, isPlotVars(i), MAPBLOCK_SCRATCH)
  end do

  !get the max and minimum variables
  if (ANY(isPlotVars)) then
     call io_getVarExtrema(NSCRATCH_GRID_VARS, globalVarMin, globalVarMax, SCRATCH)
  end if

  dowrite = 0
  if(localNumBlocks > 0) dowrite = 1
  
  !write the scratch grid vars if the user defines any in flash.par
  !we can use the same routine as when writing the unknowns.
  do i = SCRATCH_GRID_VARS_BEGIN,SCRATCH_GRID_VARS_END
     isPlotVar = isPlotVars(i)
     if(isPlotVar) then

        if(io_doublePrecision .or. io_plotfileGridQuantityDP) then
           
           unkBuf(1:NXB,1:NYB,1:NZB,1,1:localNumBlocks) = &
                scratch(i,io_ilo:io_ihi, io_jlo:io_jhi, io_klo:io_khi,1:localNumBlocks) 
           call io_h5write_unknowns(io_globalMe, &
                fileID, & 
                NXB, & 
                NYB, & 
                NZB, & 
                unkBuf, &
                globalVarMin(i), &
                globalVarMax(i), &
                io_plotGridVarStr(i), &
                localNumBlocks, &
                io_splitNumBlks,  & 
                localOffset, &
                dowrite)
           
        else
           
           unkt(1:NXB,1:NYB,1:NZB,1,1:localNumBlocks) = &
                real(scratch(i,io_ilo:io_ihi, io_jlo:io_jhi, io_klo:io_khi,1:localNumBlocks), kind = single) 
           
           spMax = real(globalVarMax(i), kind = single)
           spMin = real(globalVarMin(i), kind = single)


           call io_h5write_unknowns_sp(io_globalMe, &
                fileID, & 
                NXB,   & 
                NYB,   & 
                NZB, & 
                spMin, &
                spMax, &
                unkt,          & 
                io_plotGridVarStr(i),  & 
                localNumBlocks,  & 
                io_splitNumBlks,  & 
                localOffset, &
                dowrite)
           
        end if
     end if
  end do

  if (io_doublePrecision .or. io_plotfileGridQuantityDP) then
     deallocate(unkBuf)
  else
     deallocate(unkt)
  end if

  deallocate(isPlotVars)
  deallocate(globalVarMin)
  deallocate(globalVarMax)

#ifdef FLASH_GRID_PARAMESH3OR4
#if(NFACE_VARS>0)

    if(io_doublePrecision) then     
       
       if (io_chkGuardCellsOutput) then ! blkLimitsFX is not used otherwise
          blockID = 1
          call Grid_getBlkIndexLimits(blockID, blkLimitsFX, blkLimitsGcFX,FACEX)
#ifndef FL_NON_PERMANENT_GUARDCELLS
          blkLimitsFX = blkLimitsGcFX
#endif
       else
          allocate(faceXBuf(1,NXB+1, NYB, NZB, MAXBLOCKS))
       end if
       allocate(globalFaceXMin(NFACE_VARS))
       allocate(globalFaceXMax(NFACE_VARS))
       call io_getVarExtrema(NFACE_VARS,globalFaceXMin,globalFaceXMax,FACEX)
       
       if(NDIM .gt. 1) then
          if (io_chkGuardCellsOutput) then ! blkLimitsFY is not used otherwise
             blockID = 1
             call Grid_getBlkIndexLimits(blockID, blkLimitsFY, blkLimitsGcFY,FACEY)
#ifndef FL_NON_PERMANENT_GUARDCELLS
             blkLimitsFY = blkLimitsGcFY
#endif
          else
             allocate(faceYBuf(1,NXB, NYB+1, NZB, MAXBLOCKS))
          end if
          allocate(globalFaceYMin(NFACE_VARS))
          allocate(globalFaceYMax(NFACE_VARS))
          call io_getVarExtrema(NFACE_VARS,globalFaceYMin,globalFaceYMax,FACEY)
       end if

       if(NDIM .gt. 2) then
          if (io_chkGuardCellsOutput) then ! blkLimitsFZ is not used otherwise
             blockID = 1
             call Grid_getBlkIndexLimits(blockID, blkLimitsFZ, blkLimitsGcFZ,FACEZ)
#ifndef FL_NON_PERMANENT_GUARDCELLS
             blkLimitsFZ = blkLimitsGcFZ
#endif
          else
             allocate(faceZBuf(1,NXB, NYB, NZB+1, MAXBLOCKS))
          end if
          allocate(globalFaceZMin(NFACE_VARS))
          allocate(globalFaceZMax(NFACE_VARS))
          call io_getVarExtrema(NFACE_VARS,globalFaceZMin,globalFaceZMax,FACEZ)
       end if

       if (.NOT. io_chkGuardCellsOutput) then
       do i = 1,NFACE_VARS
          faceXBuf(1,1:NXB+1,1:NYB,1:NZB,1:MAXBLOCKS) = &
            facevarx(i,io_ilo:io_ihi+1,io_jlo:io_jhi,io_klo:io_khi,1:MAXBLOCKS)
          call io_h5write_unknowns(io_globalMe, &
               fileID, &
               NXB+1, &
               NYB, &
               NZB, &
               faceXBuf, &
               globalFaceXMin(i), &
               globalFaceXMax(i), &
               io_faceXVarLabels(i), &
               localNumBlocks, &
               io_splitNumBlks, &
               localOffset, &
               dowrite)

          if(NDIM .gt. 1) then
              
            faceYBuf(1,1:NXB,1:NYB+1,1:NZB,1:MAXBLOCKS) = &
              facevary(i,io_ilo:io_ihi,io_jlo:io_jhi+1,io_klo:io_khi,1:MAXBLOCKS)

            call io_h5write_unknowns(io_globalMe, &
                 fileID, &
                 NXB, &
                 NYB+1, &
                 NZB, &
                 faceYBuf, &
                 globalFaceYMin(i), &
                 globalFaceYMax(i), &
                 io_faceYVarLabels(i), &
                 localNumBlocks, &
                 io_splitNumBlks, &
                 localOffset, &
                 dowrite)
    
            if(NDIM .gt. 2) then
               
              faceZBuf(1,1:NXB,1:NYB,1:NZB+1,1:MAXBLOCKS) = &
               facevarz(i,io_ilo:io_ihi,io_jlo:io_jhi,io_klo:io_khi+1,1:MAXBLOCKS)
          
              call io_h5write_unknowns(io_globalMe, &
                fileID, &
                NXB, &
                NYB, &
                NZB+1, &
                faceZBuf, &
                globalFaceZMin(i), &
                globalFaceZMax(i), &
                io_faceZVarLabels(i), &
                localNumBlocks, &
                io_splitNumBlks, &
                localOffset, &
                dowrite)

            end if
    
          end if

       end do
       else ! if (.NOT. io_chkGuardCellsOutput)
          allocate(faceXBuf(1,GRID_IHI_GC+1,GRID_JHI_GC,GRID_KHI_GC,localNumBlocks))
          faceXBuf = 0.0
          do i = 1,NFACE_VARS
             faceXBuf(1,blkLimitsFX(LOW,IAXIS):blkLimitsFX(HIGH,IAXIS), &
                  blkLimitsFX(LOW,JAXIS):blkLimitsFX(HIGH,JAXIS), &
                  blkLimitsFX(LOW,KAXIS):blkLimitsFX(HIGH,KAXIS),1:localNumBlocks) = &
                  facevarx(i,:,:,:,1:localNumBlocks)
             call io_h5write_unknowns(io_globalMe, &
                  fileID, &
                  GRID_IHI_GC+1, &
                  GRID_JHI_GC, &
                  GRID_KHI_GC, &
                  faceXBuf, &
                  globalFaceXMin(i), &
                  globalFaceXMax(i), &
                  io_faceXVarLabels(i), &
                  localNumBlocks, &
                  io_splitNumBlks, &
                  localOffset, &
                  dowrite)
          end do

          if(NDIM .gt. 1) then
             allocate(faceYBuf(1,GRID_IHI_GC,GRID_JHI_GC+1,GRID_KHI_GC,localNumBlocks))
             faceYBuf = 0.0
             do i = 1,NFACE_VARS
                faceYBuf(1,blkLimitsFY(LOW,IAXIS):blkLimitsFY(HIGH,IAXIS), &
                     blkLimitsFY(LOW,JAXIS):blkLimitsFY(HIGH,JAXIS), &
                     blkLimitsFY(LOW,KAXIS):blkLimitsFY(HIGH,KAXIS),1:localNumBlocks) = &
                     facevary(i,:,:,:,1:localNumBlocks)

                call io_h5write_unknowns(io_globalMe, &
                     fileID, &
                     GRID_IHI_GC, &
                     GRID_JHI_GC+1, &
                     GRID_KHI_GC, &
                     faceYBuf, &
                     globalFaceYMin(i), &
                     globalFaceYMax(i), &
                     io_faceYVarLabels(i), &
                     localNumBlocks, &
                     io_splitNumBlks, &
                     localOffset, &
                     dowrite)
             end do

             if(NDIM .gt. 2) then

                allocate(faceZBuf(1,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC+1,localNumBlocks))
                faceZBuf = 0.0
                do i = 1,NFACE_VARS
                   faceZBuf(1,blkLimitsFZ(LOW,IAXIS):blkLimitsFZ(HIGH,IAXIS), &
                        blkLimitsFZ(LOW,JAXIS):blkLimitsFZ(HIGH,JAXIS), &
                        blkLimitsFZ(LOW,KAXIS):blkLimitsFZ(HIGH,KAXIS),1:localNumBlocks) = &
                        facevarz(i,:,:,:,1:localNumBlocks)
                   call io_h5write_unknowns(io_globalMe, &
                        fileID, &
                        GRID_IHI_GC, &
                        GRID_JHI_GC, &
                        GRID_KHI_GC+1, &
                        faceZBuf, &
                        globalFaceZMin(i), &
                        globalFaceZMax(i), &
                        io_faceZVarLabels(i), &
                        localNumBlocks, &
                        io_splitNumBlks, &
                        localOffset, &
                        dowrite)
                end do
             end if

          end if                !if (NDIM .gt. 1)
       end if                   !if (.NOT. io_chkGuardCellsOutput)

       deallocate(faceXBuf)
       deallocate(globalFaceXMin)
       deallocate(globalFaceXMax)
       
       if(NDIM .gt. 1) then
         deallocate(faceYBuf)
         deallocate(globalFaceYMin)
         deallocate(globalFaceYMax)
       end if
      
       if(NDIM .gt. 2) then 
         deallocate(faceZBuf)
         deallocate(globalFaceZMin)
         deallocate(globalFaceZMax)
       end if
       
    end if       
#endif
#endif

  if (io_globalMe .EQ. MASTER_PE) then
     write (numToStr, "(I7)") gr_globalNumBlocks
     buff = "wrote " // numToStr // " blocks"
     call Logfile_stamp( buff, '[io_writeData]')
  end if




  return

end subroutine io_writeData
