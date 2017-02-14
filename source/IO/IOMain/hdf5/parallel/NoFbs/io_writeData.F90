!!****if* source/IO/IOMain/hdf5/parallel/NoFbs/io_writeData
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

!!REORDER(5): unk, scratch, facevar[xyz]

subroutine io_writeData(fileID) 

  use IO_data, ONLY : io_globalMe, io_realParmNames, io_realParmValues, io_numRealParms, &
       io_intParmNames, io_intParmValues, io_numIntParms, &
       io_logParmNames, io_logParmValues, io_numLogParms, &
       io_strParmNames, io_strParmValues, io_numStrParms, &
       io_realScalarNames, io_realScalarValues, io_numRealScalars, &
       io_intScalarNames, io_intScalarValues, io_numIntScalars, &
       io_logScalarNames, io_logScalarValues, io_numLogScalars, &
       io_strScalarNames, io_strScalarValues, io_numStrScalars, &
       io_logToIntScalarValues, io_logToIntParmValues,  &
       io_geometry, io_setupCall, io_buildDir, io_flashRelease, &
       io_fileCreationTime, io_buildDate, io_buildMachine, io_cflags, io_fflags, &
       io_setupTimeStamp, io_buildTimeStamp, io_doublePrecision, &
       io_unklabelsGlobal, &
       io_faceXVarLabels,io_faceYVarLabels,io_faceZVarLabels,&
       io_plotVarStr, io_nPlotVars, io_bytePack, io_outputSplitNum, &
       io_ilo, io_ihi, io_jlo, io_jhi, io_klo, io_khi, io_plotGridVarStr, &
       io_iguard, io_jguard, io_kguard, io_plotfileMetadataDP, &
       io_plotfileGridQuantityDP, io_fileFormatVersion, io_meshMe, io_acrossMe, &
       io_unkNonRep, tree_data_t
  use Simulation_interface, ONLY : Simulation_mapStrToInt
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Grid_interface, ONLY : Grid_getGlobalIndexLimits, &
    Grid_getLocalNumBlks, Grid_getBlkIndexLimits, Grid_getBlkCornerID
  
  use Grid_data, ONLY : gr_gid, gr_axisMe, gr_procGrid=>gr_axisNumProcs, scratch
  
  use physicalData, only : unk,facevarx,facevary,facevarz
  use IO_interface, ONLY : IO_setScalar
  use io_typeInterface, ONLY : io_xfer_tree_data

  implicit none
  
#include "Flash_mpi.h"
#include "constants.h"
#include "Flash.h"
#include "io_flash.h"

  integer, intent(in) :: fileID
  type (tree_data_t) :: tree_data

  integer :: isize, jsize, ksize, temp
  integer :: localNumBlocks, blocksPerFile
  integer :: i, u, stat, globalNumBlocks, globalOffset
  logical :: isPlotVar

  ! allocate storage to hold a single variable information
  ! this should only be a small memory overhead
  integer, parameter :: single = SELECTED_REAL_KIND(p=6)
  !real (kind=single) :: unkt_crn(1,NXB+1,NYB+K2D,NZB+K3D,1)
  real (kind=single), allocatable :: unkt(:,:,:,:,:)
  
  ! allocate storage to hold the coordinate information and bounding box
  ! information
  real (kind=single) :: tmpSingle(MDIM,1)
  real (kind=single) :: boundBoxSingle(2,MDIM)
  
  real (kind=single) :: spMax, spMin

  real :: boundBox(2, MDIM)
  real :: blockCenterCoords(MDIM, 1) !!1 because only 1 block per proc
  real,allocatable :: unkBuf(:,:,:,:,:)
  real,allocatable :: faceXBuf(:,:,:,:,:)
  real,allocatable :: faceYBuf(:,:,:,:,:)
  real,allocatable :: faceZBuf(:,:,:,:,:)

  real, allocatable :: globalVarMin(:), globalVarMax(:)
  real, allocatable :: globalFaceXMin(:), globalFaceXMax(:)
  real, allocatable :: globalFaceYMin(:), globalFaceYMax(:)
  real, allocatable :: globalFaceZMin(:), globalFaceZMax(:)

  integer :: cornerID(MDIM), stride(MDIM)

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer, dimension(MDIM)   :: globalIndexLimits
  integer :: localnxb, localnyb, localnzb
  integer :: nxbOffset, nybOffset, nzbOffset

  integer, parameter :: presentDims = MDIM
  integer, target, dimension(1) :: nodetype = (/LEAF/), lrefine = (/1/)
  real, target, dimension(2,MDIM,1) :: bnd_box
  real, target, dimension(MDIM,1) :: coord, bsize
  
  integer :: dowrite
  
  !for nonfixed block size io we are faking that we only
  !have 1 block.  Each processor fits its portion of data
  !into that 1 block.  It also means that we have to assign
  !the value of bnd box, gid, blksize, coords as if there were
  !only one block

  !use entire physical domain to 'fake' single block
  call RuntimeParameters_get("xmin", bnd_box(LOW,IAXIS,1))
  call RuntimeParameters_get("xmax", bnd_box(HIGH,IAXIS,1))

  call RuntimeParameters_get("ymin", bnd_box(LOW,JAXIS,1))
  call RuntimeParameters_get("ymax", bnd_box(HIGH,JAXIS,1))

  call RuntimeParameters_get("zmin", bnd_box(LOW,KAXIS,1))
  call RuntimeParameters_get("zmax", bnd_box(HIGH,KAXIS,1))
  
  !get the global index limits of the domain and set the
  !values of nxb, nyb and nzb in the scalar list.
  call Grid_getGlobalIndexLimits(globalIndexLimits)

  call IO_setScalar("nxb", globalIndexLimits(IAXIS))
  call IO_setScalar("nyb", globalIndexLimits(JAXIS))
  call IO_setScalar("nzb", globalIndexLimits(KAXIS))
  call IO_setScalar("globalNumBlocks", 1)
  
  !! call the generic function prepareLists to allocate and 
  !! fill the runtime parameter lists and the scalar lists
  call io_prepareListsWrite()
  ! localNumBlocks should be 1
  call Grid_getLocalNumBlks(localNumBlocks)
    
  ! Allocate space for unkt and unkbuf
!!  allocate(unkt(1,NXB,NYB,NZB,1))
 
 !get size of block, this could change from block to block so we can't use NXB, NYB, NZB 
  call Grid_getBlkIndexLimits(1, blkLimits, blkLimitsGC, CENTER)
  localnxb = blkLimits(HIGH,IAXIS) - blkLimits(LOW,IAXIS) + 1
  localnyb = blkLimits(HIGH,JAXIS) - blkLimits(LOW,JAXIS) + 1
  localnzb = blkLimits(HIGH,KAXIS) - blkLimits(LOW,KAXIS) + 1
  
  if (io_doublePrecision .or. io_plotfileGridQuantityDP) then 
    allocate(unkBuf(1,localnxb,localnyb,localnzb,1))
  else 
    allocate(unkt(1, localnxb, localnyb, localnzb, 1))
  end if

  call Grid_getBlkCornerID(1, cornerID, stride)

  nxbOffset = cornerID(IAXIS) -1
  nybOffset = (cornerID(JAXIS) -1) * K2D
  nzbOffset = (cornerID(KAXIS) -1) * K3D
  
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
  


  !for this non fixed block size io only the master processor will write out the 
  !tree data
  if(io_meshMe == MASTER_PE .and. io_acrossMe == 0) then
     localNumBlocks=1
  else
     localNumBlocks = 0
  end if
  globalNumBlocks=1
  globalOffset=0

  call io_createDatasets(fileID, globalNumBlocks, presentDims)

  !if faking one block, everything is a boundary, no neighbors
  gr_gid = -21
  !use entire physical domain to 'fake' single block
  coord(:,1) = bnd_box(HIGH,:,1) * 0.5
  bsize(:,1) = bnd_box(HIGH,:,1) - bnd_box(LOW,:,1)

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




  !get maximum and minimum for all variables
  allocate(globalVarMin(ubound(io_unklabelsGlobal,1)))
  allocate(globalVarMax(ubound(io_unklabelsGlobal,1)))
  
#if(NFACE_VARS > 0)
  allocate(globalFaceXMin(NFACE_VARS))
  allocate(globalFaceXMax(NFACE_VARS))
  if(NDIM > 1) then
    allocate(globalFaceYMin(NFACE_VARS))
    allocate(globalFaceYMax(NFACE_VARS))
  end if !NDIM > 1
  if(NDIM > 2) then
    allocate(globalFaceZMin(NFACE_VARS))
    allocate(globalFaceZMax(NFACE_VARS))
  end if !NDIM > 2
#endif
  
  call io_getVarExtrema(ubound(io_unklabelsGlobal,1), globalVarMin, globalVarMax, CENTER)

#if(NFACE_VARS > 0)
  call io_getVarExtrema(NFACE_VARS,globalFaceXMin,globalFaceXMax, FACEX)
  if (NDIM > 1) call io_getVarExtrema(NFACE_VARS,globalFaceYMin,globalFaceYMax, FACEY)
  if (NDIM > 2) call io_getVarExtrema(NFACE_VARS,globalFaceZMin,globalFaceZMax, FACEZ)
#endif
  !--------------------------------------------------------------------------
  ! store the unknowns -- 
  !--------------------------------------------------------------------------
!  do i = UNK_VARS_BEGIN,UNK_VARS_END
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
       unkBuf(1,1:localnxb,1:localnyb,1:localnzb,1) = &           
            unk(i,1+io_iguard:localnxb+io_iguard, 1+io_jguard:localnyb+io_jguard, 1+io_kguard:localnzb+io_kguard,1)
      call io_h5write_unknowns(io_globalMe, &
        fileID, &
        globalIndexLimits(IAXIS), &
        globalIndexLimits(JAXIS), &
        globalIndexLimits(KAXIS), &
        nxbOffset, &
        nybOffset, &
        nzbOffset, &
        localnxb, & 
        localnyb, & 
        localnzb, & 
        unkBuf, & 
        globalVarMin(u), &
        globalVarMax(u), &
        io_unklabelsGlobal(u), &
        dowrite)
    else
      !FUTURE: still need to handle corner data at some point
      isPlotVar = any(io_unklabelsGlobal(u) == io_plotVarStr(1:io_nPlotVars))
      if (isPlotVar) then
         if (io_plotfileGridQuantityDP) then
            unkBuf(1,1:localnxb,1:localnyb,1:localnzb,1) = &           
                 unk(i,1+io_iguard:localnxb+io_iguard, 1+io_jguard:localnyb+io_jguard, 1+io_kguard:localnzb+io_kguard,1)
            call io_h5write_unknowns(io_globalMe, &
                 fileID, &
                 globalIndexLimits(IAXIS), &
                 globalIndexLimits(JAXIS), &
                 globalIndexLimits(KAXIS), &
                 nxbOffset, &
                 nybOffset, &
                 nzbOffset, &
                 localnxb, & 
                 localnyb, & 
                 localnzb, & 
                 unkBuf, & 
                 globalVarMin(u), &
                 globalVarMax(u), &
                 io_unklabelsGlobal(u), &
                 dowrite)
         else
            unkt(1,:,:,:,1) = real(unk(i,1+io_iguard:localnxb+io_iguard, &
                 1+io_jguard:localnyb+io_jguard, &
                 1+io_kguard:localnzb+io_kguard, &
                 1), kind = single)
            call io_h5write_unknowns_sp(io_globalMe, &
                 fileID, & 
                 globalIndexLimits(IAXIS), &
                 globalIndexLimits(JAXIS), &
                 globalIndexLimits(KAXIS), &
                 nxbOffset, &
                 nybOffset, &
                 nzbOffset, &
                 localnxb,   & 
                 localnyb,   & 
                 localnzb, & 
                 unkt, & 
                 real (globalVarMin(u), kind=single), &
                 real (globalVarMax(u), kind=single), &
                 io_unklabelsGlobal(u), &
                 dowrite)
         end if
      end if
    end if
  end do

  if(io_doublePrecision .or. io_plotfileGridQuantityDP) then
    deallocate(unkBuf)
  else
    deallocate(unkt)
  end if

  !!Write out facevars
  !!DEV: Face-centered variables are not sent to plot files.
#if(NFACE_VARS > 0) 
  allocate(faceXBuf(1,localnxb+1,localnyb,localnzb,1))
  if(NDIM > 1) allocate(faceYBuf(1,localnxb,localnyb+1,localnzb,1))
  if(NDIM > 2) allocate(faceZBuf(1,localnxb,localnyb,localnzb+1,1))

  ! all facevars are replicated
  dowrite = 0
  if(io_acrossMe .eq. 0) dowrite = 1
  
  do i = 1, NFACE_VARS
    faceXBuf(1,1:localnxb+1,1:localnyb,1:localnzb,1) = &
      facevarx(i, &
               1+io_iguard:localnxb+io_iguard+1, &
               1+io_jguard:localnyb+io_jguard, &
               1+io_kguard:localnzb+io_kguard, &
               1)
    
    if(gr_axisMe(IAXIS) == (gr_procGrid(IAXIS)- 1)) then
      call io_h5write_unknowns(io_globalMe,&
                               fileID, &
                               globalIndexLimits(IAXIS) + 1, &
                               globalIndexLimits(JAXIS), &
                               globalIndexLimits(KAXIS), &
                               nxbOffset, & 
                               nybOffset, &
                               nzbOffset, &
                               localnxb + 1, &
                               localnyb, &
                               localnzb, &
                               faceXBuf(1,1:localnxb+1,:,:,1), &
                               globalFaceXMin(i), &
                               globalFaceXMax(i), &
                               io_faceXVarLabels(i), &
                               dowrite)
    else
      call io_h5write_unknowns(io_globalMe, &
                               fileID, &
                               globalIndexLimits(IAXIS) + 1, &
                               globalIndexLimits(JAXIS), &
                               globalIndexLimits(KAXIS), &
                               nxbOffset, & 
                               nybOffset, &
                               nzbOffset, &
                               localnxb, &
                               localnyb, &
                               localnzb, &
                               faceXBuf(1,1:localnxb,:,:,1), &
                               globalFaceXMin(i), &
                               globalFaceXMax(i), &
                               io_faceXVarLabels(i), &
                               dowrite)
    end if

    if(NDIM > 1) then
!!$
      faceYBuf(1,1:localnxb,1:localnyb+1,1:localnzb,1) = &
        facevary(i, &
                 1+io_iguard:localnxb+io_iguard, &
                 1+io_jguard:localnyb+io_jguard+1, &
                 1+io_kguard:localnzb+io_kguard, &
                 1)
      if(gr_axisMe(JAXIS) == (gr_procGrid(JAXIS)-1)) then
        call io_h5write_unknowns(io_globalMe, &
                                 fileID, &
                                 globalIndexLimits(IAXIS), &
                                 globalIndexLimits(JAXIS) + 1, &
                                 globalIndexLimits(KAXIS), &
                                 nxbOffset, & 
                                 nybOffset, &
                                 nzbOffset, &
                                 localnxb, &
                                 localnyb + 1, &
                                 localnzb, &
                                 faceYBuf(1,:,1:localnyb+1,:,1), &
                                 globalFaceYMin(i), &
                                 globalFaceYMax(i), &
                                 io_faceYVarLabels(i), &
                                 dowrite)
      else
        call io_h5write_unknowns(io_globalMe, &
                                 fileID, &
                                 globalIndexLimits(IAXIS), &
                                 globalIndexLimits(JAXIS) + 1, &
                                 globalIndexLimits(KAXIS), &
                                 nxbOffset, & 
                                 nybOffset, &
                                 nzbOffset, &
                                 localnxb, &
                                 localnyb, &
                                 localnzb, &
                                 faceYBuf(1,:,1:localnyb,:,1), &
                                 globalFaceYMin(i), &
                                 globalFaceYMax(i), &
                                 io_faceYVarLabels(i), &
                                 dowrite)
      end if
    end if !NDIM > 1

    if(NDIM > 2) then
      faceZBuf(1,1:localnxb,1:localnyb,1:localnzb+1,1) = &
        facevarz(i, &
                 1+io_iguard:localnxb+io_iguard, &
                 1+io_jguard:localnyb+io_jguard, &
                 1+io_kguard:localnzb+io_kguard+1, &
                 1)
      if(gr_axisMe(KAXIS) == (gr_procGrid(KAXIS)-1)) then
        call io_h5write_unknowns(io_globalMe, &
                                 fileID, &
                                 globalIndexLimits(IAXIS), &
                                 globalIndexLimits(JAXIS), &
                                 globalIndexLimits(KAXIS) + 1, &
                                 nxbOffset, & 
                                 nybOffset, &
                                 nzbOffset, &
                                 localnxb, &
                                 localnyb, &
                                 localnzb + 1, &
                                 faceZBuf(1,:,:,1:localnzb+1,1), &
                                 globalFaceZMin(i), &
                                 globalFaceZMax(i), &
                                 io_faceZVarLabels(i), &
                                 dowrite)
      else
        call io_h5write_unknowns(io_globalMe, &
                                 fileID, &
                                 globalIndexLimits(IAXIS), &
                                 globalIndexLimits(JAXIS), &
                                 globalIndexLimits(KAXIS) + 1, &
                                 nxbOffset, & 
                                 nybOffset, &
                                 nzbOffset, &
                                 localnxb, &
                                 localnyb, &
                                 localnzb, &
                                 faceZBuf(1,:,:,1:localnzb,1), &
                                 globalFaceZMin(i), &
                                 globalFaceZMax(i), &
                                 io_faceZVarLabels(i), &
                                 dowrite)
      end if
    end if !NDIM > 2
  end do
  
  deallocate(faceXBuf)
  if(NDIM > 1) deallocate(faceYBuf)
  if(NDIM > 2) deallocate(faceZBuf)
#endif
  
  deallocate(globalVarMin)
  deallocate(globalVarMax)
#if(NFACE_VARS > 0)
  deallocate(globalFaceXMin)
  deallocate(globalFaceXMax)
  if(NDIM > 1) then
    deallocate(globalFaceYMin)
    deallocate(globalFaceYMax)
  end if
  if(NDIM > 2) then
    deallocate(globalFaceZMin)
    deallocate(globalFaceZMax)
  end if
#endif

 end subroutine io_writeData


!  allocate(unkBuf(NXB,NYB,NZB,1,1))
!  allocate(unkt(NXB,NYB,NZB,1,1))
!
!  !write the scratch grid vars if the user defines any in flash.par
!  !we can use the same routine as when writing the unknowns.
!  do i = SCRATCH_GRID_VARS_BEGIN,SCRATCH_GRID_VARS_END
!
!     call io_isPlotVar(i, isPlotVar, MAPBLOCK_SCRATCH)
!     if(isPlotVar) then
!        
!        if(io_doublePrecision) then
!           
!           unkBuf(1:NXB,1:NYB,1:NZB,1,1) = scratch(io_ilo:io_ihi, io_jlo:io_jhi, io_klo:io_khi,i,1) 
!           call io_h5write_unknowns(io_globalMe, &
!                fileID, & 
!                NXB, & 
!                NYB, & 
!                NZB, & 
!                unkBuf, & 
!                io_plotGridVarStr(i), &
!                localNumBlocks, &
!                globalNumBlocks,  & 
!                globalOffset)
!           
!        else
!           
!           unkt(1:NXB,1:NYB,1:NZB,1,1) = real(scratch(io_ilo:io_ihi, io_jlo:io_jhi, io_klo:io_khi,i,1), kind = single) 
!           
!           call io_h5write_unknowns_sp(io_globalMe, &
!                fileID, & 
!                NXB,   & 
!                NYB,   & 
!                NZB, & 
!                spMax, &
!                spMin, &
!                unkt,          & 
!                io_plotGridVarStr(i),  & 
!                localNumBlocks,  & 
!                globalNumBlocks,  & 
!                globalOffset)
!           
!        end if
!     end if
!  end do
!
!  deallocate(unkBuf)
!  deallocate(unkt)
!
!
