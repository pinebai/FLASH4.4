!!****if* source/IO/IOMain/hdf5/parallel/UG/io_writeData
!!
!! NAME
!!
!!  io_writeData
!!
!!
!! SYNOPSIS
!!
!!  io_writeData(int(in) :: fileID) 
!!                
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
!!  A single record for each of the data structures is created.  A
!!  processor only writes to a subset of this record.  Each record has a
!!  dimension with length = tot_blocks.  The offset of a processor into this
!!  dimension is computed by looking at the total number of blocks that are
!!  below the current processor.
!!
!!  In this version of io_writeData, each variable is given its own
!!  record -- this makes it easier to change the variable list in the
!!  future without disturbing the format of the file.
!!
!!  The include file -- hdf5_flash.h is used for the C routines 
!!
!!  
!!
!! ARGUMENTS
!! 
!!  fileID - integer file identifier for hdf5 file
!!
!! NOTES
!!
!!  This current version of the Uniform Grid IO "fakes" some paramesh or
!!  tree structure data.  For example, the uniform grid is all at the same
!!  resolution so the concept of refinement levels doesn't really exist.  
!!  However, for IO purposes and to make our visualization packages uniform
!!  all checkpoint files must write out a refinement level.  For the Uniform
!!  Grid, this is simply set to 1.  The same thing goes for node type, Paramesh
!!  and other grid packages need to keep track of whether a block is a parent
!!  leaf or ancestor.  The Uniform Grid sets all blocks to type LEAF.  This 
!!  is also necessary in case the a user wants to perform some initial 
!!  conditions with the Uniform Grid and then perform the bulk of the
!!  calculation with an amr package.
!!
!!
!!***

!!REORDER(5):unk, unkBuf, facevar[xyz], face[XYZ]Buf


subroutine io_writeData(fileID) 

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
       io_unklabelsGlobal, io_plotVarStr, io_nPlotVars, io_outputSplitNum, &
       io_ilo, io_ihi, io_jlo, io_jhi, io_klo, io_khi, io_plotGridVarStr,&
       io_faceXVarLabels, io_faceYVarLabels, io_faceZVarLabels, &
       io_plotfileMetadataDP, io_plotfileGridQuantityDP, io_fileFormatVersion, &
       io_meshMe, io_meshNumProcs, io_acrossMe, io_unkNonRep, tree_data_t
  use Simulation_interface, ONLY : Simulation_mapStrToInt
  
  use Grid_interface, ONLY : Grid_getLocalNumBlks, &
    Grid_getBlkBoundBox, Grid_getBlkCenterCoords, &
    Grid_getBlkPhysicalSize

  use Grid_data, ONLY : gr_gid,scratch,scratch_ctr,scratch_facevarx,&
       scratch_facevary,scratch_facevarz
  
  use physicalData, only : unk,  facevarx, facevary, facevarz

  use io_typeInterface, ONLY : io_xfer_tree_data
  implicit none
  
#include "Flash_mpi.h"
#include "constants.h"
#include "Flash.h"
#include "io_flash.h"

  integer, intent(in) ::  fileID

  type (tree_data_t) :: tree_data
  integer :: isize, jsize, ksize
  integer :: localNumBlocks
  integer :: i, u, stat, globalNumBlocks, globalOffset
  logical :: isPlotVar

  ! allocate storage to hold a single variable information
  ! this should only be a small memory overhead
  integer, parameter :: single = SELECTED_REAL_KIND(p=6)
  !real (kind=single) :: unkt_crn(1,NXB+1,NYB+K2D,NZB+K3D,1)
  real (kind=single), allocatable :: unkt(:,:,:,:,:)
  integer :: blkLimits(HIGH,MDIM), blkLimitsGC(HIGH,MDIM)
  
  ! allocate storage to hold the coordinate information and bounding box
  ! information
  real (kind=single) :: tmpSingle(MDIM)
  real (kind=single) :: boundBoxSingle(2,MDIM)
  
  real (kind=single) :: spMax, spMin
  real :: blockSize(MDIM)
  real :: boundBox(2, MDIM)
  real :: blockCenterCoords(MDIM)
  real,allocatable :: unkBuf(:,:,:,:,:)
  real,allocatable :: faceXBuf(:,:,:,:,:)
  real,allocatable :: faceYBuf(:,:,:,:,:)
  real,allocatable :: faceZBuf(:,:,:,:,:)

  logical, allocatable :: isPlotVars(:)
  real, allocatable :: globalVarMin(:), globalVarMax(:)
  real, allocatable :: faceXVarMin(:), faceXVarMax(:)
  real, allocatable :: faceYVarMin(:), faceYVarMax(:)
  real, allocatable :: faceZVarMin(:), faceZVarMax(:)

  integer :: dowrite

  integer, target, dimension(1) :: nodetype = (/LEAF/), lrefine = (/1/)
  real, target, dimension(2,MDIM,1) :: bnd_box
  real, target, dimension(MDIM,1) :: coord, bsize
  integer, parameter :: presentDims = MDIM

  !! call the generic function prepareLists to allocate and 
  !! fill the runtime parameter lists and the scalar lists
  call io_prepareListsWrite()
    
  ! localNumBlocks should be 1
  call Grid_getLocalNumBlks(localNumBlocks)
  
  ! Allocate space for unkt and unkbuf

  allocate(unkBuf(1,NXB,NYB,NZB,1))
  allocate(unkt(1,NXB,NYB,NZB,1))


  globalNumBlocks = io_meshNumProcs !io_globalNumProcs
  
  globalOffset = io_meshMe !io_globalMe


  
  if(io_doublePrecision) then

     call io_h5write_header(io_meshMe, ubound(io_unklabelsGlobal,1), fileID, io_geometry, &
          io_unklabelsGlobal,io_setupCall, io_fileCreationTime, io_flashRelease, &
          io_buildDate, io_buildDir, io_buildMachine, io_cflags, io_fflags, &
          io_setupTimeStamp, io_buildTimeStamp, io_fileFormatVersion, &
          io_outputSplitNum)

  else

     call io_h5write_header(io_meshMe, io_nPlotVars, fileID, io_geometry, &
          io_plotVarStr,io_setupCall, io_fileCreationTime, io_flashRelease, &
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



  call io_createDatasets(fileID, globalNumBlocks, presentDims)

  call Grid_getBlkBoundBox(1, boundBox)
  bnd_box(:,:,1) = boundBox(:,:)
  call Grid_getBlkCenterCoords(1, blockCenterCoords)
  coord(:,1) = blockCenterCoords(:)
  call Grid_getBlkPhysicalSize(1, blockSize)
  bsize(:,1) = blockSize(:)

  tree_data % bnd_box => bnd_box
  tree_data % coord => coord
  tree_data % bsize => bsize
  tree_data % gid => gr_gid
  tree_data % nodetype => nodetype
  tree_data % lrefine => lrefine
  nullify(tree_data % bflags)
  nullify(tree_data % which_child)
  allocate(tree_data % procnumber(max(1,localNumBlocks)))
  tree_data % procnumber(:) = io_meshMe
  nullify(tree_data % gsurr_blks)

  call io_xfer_tree_data(tree_data, fileID, IO_FILE_HDF5, IO_WRITE_XFER, &
       localNumBlocks, globalOffset, presentDims)

  deallocate(tree_data % procnumber)
  nullify(tree_data % procnumber)



  allocate(globalVarMin(ubound(io_unklabelsGlobal,1)))
  allocate(globalVarMax(ubound(io_unklabelsGlobal,1)))

  !get the max and minimum variables
  call io_getVarExtrema(ubound(io_unklabelsGlobal,1), globalVarMin, globalVarMax, CENTER)

#if (NFACE_VARS > 0)
      allocate(faceXVarMin(NFACE_VARS))
      allocate(faceXVarMax(NFACE_VARS))
      call io_getVarExtrema(NFACE_VARS, faceXVarMin, faceXVarMax, FACEX)
      if(NDIM .gt. 1) then
        allocate(faceYVarMin(NFACE_VARS))
        allocate(faceYVarMax(NFACE_VARS))
        call io_getVarExtrema(NFACE_VARS, faceYVarMin, faceYVarMax, FACEY)
      end if
      if(NDIM .gt. 2) then
        allocate(faceZVarMin(NFACE_VARS))
        allocate(faceZVarMax(NFACE_VARS))
        call io_getVarExtrema(NFACE_VARS, faceZVarMin, faceZVarMax, FACEZ)
      end if
#endif


  !--------------------------------------------------------------------------
  ! store the unknowns -- 
  !--------------------------------------------------------------------------

  !do i = UNK_VARS_BEGIN,UNK_VARS_END
  do u=1, ubound(io_unklabelsGlobal,1)
     call Simulation_mapStrToInt(io_unklabelsGlobal(u),i,MAPBLOCK_UNK)
     ! only write mesh replicated data from mesh 0
     dowrite = 0
     if(i /= NONEXISTENT) then
        if(io_acrossMe .eq. 0 .or. io_unkNonRep(i) > 0) dowrite = 1
     else
        i = 1 ! gotta give it something even though the data wont be written
     end if

     if(io_doublePrecision) then
        unkBuf(1,:,:,:,1) = unk(i,io_ilo:io_ihi, io_jlo:io_jhi, io_klo:io_khi,1) 

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
             globalNumBlocks,  & 
             globalOffset, &
             dowrite)

     else
        !check to see if variable should be written to the plotfile
        isPlotVar = any(io_unklabelsGlobal(u) == io_plotVarStr(1:io_nPlotVars))
        if (isPlotVar) then
           
           if (io_plotfileGridQuantityDP) then
              unkBuf(1,:,:,:,1) = unk(i,io_ilo:io_ihi, io_jlo:io_jhi, io_klo:io_khi,1) 
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
                   globalNumBlocks,  & 
                   globalOffset, &
                   dowrite)
           else
              unkt(1,:,:,:,1) = real(unk(i,io_ilo:io_ihi, &
                   io_jlo:io_jhi, &
                   io_klo:io_khi, &
                   1), kind = single)
              spMin = real(globalVarMin(u), kind = single)
              spMax = real(globalVarMax(u), kind = single)
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
                   globalNumBlocks,  & 
                   globalOffset, &
                   dowrite)
           end if
        end if
     end if
       
  end do

  
  deallocate(unkBuf)
  deallocate(unkt)
  allocate(unkBuf(NXB,NYB,NZB,1,1))
  allocate(unkt(NXB,NYB,NZB,1,1))


  deallocate(globalVarMin)
  deallocate(globalVarMax)
  
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

  !write the scratch grid vars if the user defines any in flash.par
  !we can use the same routine as when writing the unknowns.
  do i = SCRATCH_GRID_VARS_BEGIN,SCRATCH_GRID_VARS_END
     dowrite = 1
     isPlotVar = isPlotVars(i)
     if(isPlotVar) then
        
        if(io_doublePrecision .or. io_plotfileGridQuantityDP) then
           
           unkBuf(1:NXB,1:NYB,1:NZB,1,1) = scratch(i,io_ilo:io_ihi, io_jlo:io_jhi, io_klo:io_khi,1) 
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
                globalNumBlocks,  & 
                globalOffset, &
                dowrite)
           
        else
           
           unkt(1:NXB,1:NYB,1:NZB,1,1) = real(scratch(i,io_ilo:io_ihi, io_jlo:io_jhi, io_klo:io_khi,1), kind = single) 
           
           spMin = real(globalVarMin(i), kind = single)
           spMax = real(globalVarMax(i), kind = single)


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
                globalNumBlocks,  & 
                globalOffset, &
                dowrite)
           
        end if
     end if
  end do

  deallocate(unkBuf)
  deallocate(unkt)

  deallocate(isPlotVars)
  deallocate(globalVarMin)
  deallocate(globalVarMax)

#if (NFACE_VARS>0)

    if (io_doublePrecision) then
      allocate(faceXBuf(1,NXB+1,NYB,NZB,1))
      if(NDIM .gt. 1) allocate(faceYBuf(1,NXB,NYB+1,NZB,1))
      if(NDIM .gt. 2) allocate(faceZBuf(1, NXB,NYB,NZB+1,1))
      dowrite = 1
      
      do i = 1,NFACE_VARS
        faceXBuf(1,:,:,:,1) = facevarx(i,io_ilo:io_ihi+1, io_jlo:io_jhi, io_klo:io_khi,1) 
        call io_h5write_unknowns(io_globalMe, &
                                 fileID, &
                                 NXB + 1, &
                                 NYB, &
                                 NZB, &
                                 faceXBuf(1,:,:,:,1), &
                                 faceXVarMin(i), &
                                 faceXVarMax(i), &
                                 io_faceXVarLabels(i), &
                                 localNumBlocks, &
                                 globalNumBlocks, &
                                 globalOffset, &
                                 dowrite)
        
        if(NDIM .gt. 1) then
        
          faceYBuf(1,:,:,:,1) = facevary(i,io_ilo:io_ihi, io_jlo:io_jhi+1, io_klo:io_khi,1) 
          call io_h5write_unknowns(io_globalMe, &
                                 fileID, &
                                 NXB, &
                                 NYB + 1, &
                                 NZB, &
                                 faceYBuf(1,:,:,:,1), &
                                 faceYVarMin(i), &
                                 faceYVarMax(i), &
                                 io_faceYVarLabels(i), &
                                 localNumBlocks, &
                                 globalNumBlocks, &
                                 globalOffset, &
                                 dowrite)
        
        end if

        if(NDIM .gt. 2) then

          faceZBuf(1,:,:,:,1) = facevarz(i,io_ilo:io_ihi, io_jlo:io_jhi, io_klo:io_khi+1,1) 
          call io_h5write_unknowns(io_globalMe, &
                                   fileID, &
                                   NXB, &
                                   NYB, &
                                   NZB + 1, &
                                   faceZBuf(1,:,:,:,1), &
                                   faceZVarMin(i), &
                                   faceZVarMax(i), &
                                   io_faceZVarLabels(i), &
                                   localNumBlocks, &
                                   globalNumBlocks, &
                                   globalOffset, &
                                   dowrite)
        end if

      end do
      deallocate(faceXBuf)
      deallocate(faceXVarMin)
      deallocate(faceXVarMax)
      if(NDIM .gt. 1) then
        deallocate(faceYBuf)
        deallocate(faceYVarMin)
        deallocate(faceYVarMax)
      end if
      if(NDIM .gt. 2) then
        deallocate(faceZBuf)
        deallocate(faceZVarMin)
        deallocate(faceZVarMax)
      end if
    end if
#endif


end subroutine io_writeData
