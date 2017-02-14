!!****if* source/IO/IOMain/hdf5/serial/UG/io_writeData
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
!!  paramesh data.  IO is done in serial -- data is copied to processor 0 and then
!!  written out to file
!!  to do the writing is performed.  HDF5 v. 1.4.0 or later is required
!!
!!  HDF5 uses MPI-IO (via ROMIO) to support parallel IO.  Each processor
!!  must open the file, create the datasets, and dataspaces for each HDF
!!  record.
!!
!!  A single record for each of the PARAMESH data structures is created.  A
!!  processor only writes to a subset of this record.  Each record has a
!!  dimension with length = tot_blocks.  The offset of a processor into this
!!  dimension is computed by looking at the total number of blocks that are
!!  below the current processor.
!!
!!  In this version of io_writeData, each variable is given its own
!!  record -- this makes it easier to change the variable list in the
!!  future without disturbing the format of the file.
!!
!!  This routine is used to write both checkpoint files and plotfiles.  
!!  Checkpoint files are written in double precision, while plotfiles 
!!  are written in single precision unless io_plotfileGridQuantityDP
!!  is TRUE.
!!
!! ARGUMENTS
!! 
!!  fileID - integer file identifier for hdf5 file
!!
!!***

!!REORDER(5): unk, facevar[xyz]

subroutine io_writeData (fileID)

 use IO_data, ONLY : io_globalMe, io_globalNumProcs, io_globalComm,&
       io_realParmNames, io_realParmValues, io_numRealParms, &
       io_intParmNames, io_intParmValues, io_numIntParms, &
       io_logParmNames, io_logParmValues, io_numLogParms, &
       io_strParmNames, io_strParmValues, io_numStrParms, &
       io_realScalarNames, io_realScalarValues, io_numRealScalars, &
       io_intScalarNames, io_intScalarValues, io_numIntScalars, &
       io_logScalarNames, io_logScalarValues, io_numLogScalars, &
       io_strScalarNames, io_strScalarValues, io_numStrScalars, &
       io_logToIntScalarValues, io_logToIntParmValues, io_unklabels, &
       io_ilo, io_ihi, io_jlo, io_jhi, io_klo, io_khi, io_geometry, &
       io_setupCall, io_buildDir, io_flashRelease, &
       io_fileCreationTime, io_buildDate, io_buildMachine, io_cflags, io_fflags, &
       io_setupTimeStamp, io_buildTimeStamp, io_doublePrecision, &
       io_outputSplitNum, io_plotVarStr, io_nPlotVars, io_plotGridVarStr, &
       io_faceXVarLabels, io_faceYVarLabels, io_faceZVarLabels, &
       io_plotfileMetadataDP, io_plotfileGridQuantityDP, io_fileFormatVersion, &
       tree_data_t
 use Driver_interface, ONLY : Driver_abortFlash
 use Grid_interface, ONLY : Grid_getLocalNumBlks, &
   Grid_getBlkBoundBox, Grid_getBlkCenterCoords, &
   Grid_getBlkPhysicalSize
  
  
  use Grid_data, ONLY : gr_gid, scratch

  use physicaldata, ONLY :unk,facevarx,facevary,facevarz
  
  use io_typeInterface, ONLY : io_xfer_tree_data

  implicit none


#include "Flash_mpi.h"
#include "constants.h"
#include "Flash.h"
#include "io_flash.h"

  integer, INTENT(in) :: fileID

  type (tree_data_t) :: tree_data
  integer :: blockID = 1
  integer :: jproc, i, j
  integer :: globalNumBlocks, globalOffset, ngid 
  integer :: ierr


  ! Define the temp variables (t) which hold the data transfered
  ! to the MASTER_PE
  integer,target,allocatable :: gidt(:,:)
  integer :: status(MPI_STATUS_SIZE)

  real :: boundBox(2,MDIM)
  real :: blockSize(MDIM)
  real :: blockCenterCoords(MDIM)

  real, target :: boundBoxt(2,MDIM,1)
  real, target :: bsize(MDIM,1)
  real, target :: blockCenterCoordst(MDIM,1)

  real, allocatable :: unkt(:,:,:,:)
  real, allocatable :: scratcht(:,:,:,:)
  real, allocatable :: facext(:,:,:,:)
  real, allocatable :: faceyt(:,:,:,:)
  real, allocatable :: facezt(:,:,:,:)


  real  :: globalVarMin(NUNK_VARS), globalVarMax(NUNK_VARS)
  real  :: globalVarMinScratch(NSCRATCH_GRID_VARS), globalVarMaxScratch(NSCRATCH_GRID_VARS)
  real  :: globalFaceXMin(NUNK_VARS), globalFaceXMax(NUNK_VARS)
  real  :: globalFaceYMin(NUNK_VARS), globalFaceYMax(NUNK_VARS)
  real  :: globalFaceZMin(NUNK_VARS), globalFaceZMax(NUNK_VARS)
  
  character (len=MAX_STRING_LENGTH) :: filename

  integer :: localNumBlocks

  integer :: offset

  ! for logfile output
  character(len=MAX_STRING_LENGTH), save, allocatable, dimension(:,:) :: strBuff
  character(len=16)                                    :: num_to_str


  ! block data message buffering stuff. not using this right now
  integer, parameter :: MAX_TRANS_SIZE = 4000000  ! Maximum number of bytes
  ! to transfer at one time
  ! when sending block data
  integer, save      :: MAX_BLK_MSGS
  integer            :: n_blk_msgs, pblkcount


  ! allocate storage to hold a single variable information
  ! this should only be a small memory overhead
  integer, parameter :: single = SELECTED_REAL_KIND(p=6)

  real (kind=single) :: singleUnk(NXB,NYB,NZB,1)
  real (kind=single) :: singleScratch(NXB,NYB,NZB,MAXBLOCKS)


  ! allocate storage to hold the coordinate information and bounding box
  ! information
  real (kind=single) :: sizeSingle(MDIM)
  real (kind=single) :: coordSingle(MDIM)
  real (kind=single) :: bndSingle(2,MDIM)
  real (kind=single) :: spMax, spMin

  logical :: isPlotVar
  integer, parameter :: presentDims = MDIM
  integer, target, dimension(1) :: nodetype = (/LEAF/), lrefine = (/1/)
  integer, target, dimension(1) :: procnumber

!=============================================================================
  

  

  call Grid_getLocalNumBlks(localNumBlocks)
  if (localNumBlocks /= 1) then
     call Driver_abortFlash("Error: UG must have 1 block per proc")
  end if


  globalNumBlocks = io_globalNumProcs
  globalOffset = io_globalMe

  call Grid_getBlkBoundBox(1, boundBox)

  call Grid_getBlkCenterCoords(1, blockCenterCoords)

  call Grid_getBlkPhysicalSize(1, blockSize)
  bsize(:,1) = blockSize(:)

  if (NDIM == 1) then
        ngid = 5
     else if (NDIM == 2) then
        ngid = 9
     else if (NDIM == 3) then
        ngid = 15
  end if

  allocate(gidt(ngid, 1))


  tree_data % bnd_box => boundBoxt
  tree_data % coord => blockCenterCoordst
  tree_data % bsize => bsize
  tree_data % gid => gidt
  tree_data % nodetype => nodetype
  tree_data % lrefine => lrefine
  nullify(tree_data % bflags)
  nullify(tree_data % which_child)
  tree_data % procnumber => procnumber
  nullify(tree_data % gsurr_blks)

  
  if (io_globalMe == MASTER_PE) then
     
     if(io_doublePrecision) then
        call io_h5write_header(io_globalMe, NUNK_VARS, fileID, io_geometry, &
             io_unklabels, io_setupCall, io_fileCreationTime, io_flashRelease, &
             io_buildDate, io_buildDir, io_buildMachine, io_cflags, io_fflags, &
             io_setupTimeStamp, io_buildTimeStamp, io_fileFormatVersion, &
             io_outputSplitNum)

     else

        call io_h5write_header(io_globalMe, io_nPlotVars, fileID, io_geometry, &
             io_plotVarStr, io_setupCall, io_fileCreationTime, io_flashRelease, &
             io_buildDate, io_buildDir, io_buildMachine, io_cflags, io_fflags, &
             io_setupTimeStamp, io_buildTimeStamp, io_fileFormatVersion, &
             io_outputSplitNum)
     end if


     call io_prepareListsWrite()

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

     ! keep track of the global offset. With UG this just becomes proc number
     offset = 0

     call io_createDatasets(fileID, globalNumBlocks, presentDims)
  end if



  !get the max and minimum variables
  call io_getVarExtrema(NUNK_VARS, globalVarMin, globalVarMax, CENTER)
  call io_getVarExtrema(NSCRATCH_GRID_VARS, globalVarMinScratch, globalVarMaxScratch, SCRATCH)

!!GET FACEVAR EXTREMA
#if (NFACE_VARS > 0)
    call io_getVarExtrema(NFACE_VARS, globalFaceXMin, globalFaceXMax, FACEX)
    if(NDIM .gt. 1) call io_getVarExtrema(NFACE_VARS, globalFaceYMin, globalFaceYMax, FACEY)
    if(NDIM .gt. 2) call io_getVarExtrema(NFACE_VARS, globalFaceZMin, globalFaceZMax, FACEZ)
#endif

  !-----------------------------------------------------------------------------
  ! loop over all of the processors.  All the data is moved to processor 0 for
  ! storage using MPI sends and receives.
  !-----------------------------------------------------------------------------

  do jproc = 0,io_globalNumProcs-1

     if (io_globalMe == MASTER_PE) then

        if (jproc /= 0) then
           
           call MPI_RECV(blockCenterCoordst(1,1), MDIM, FLASH_REAL, &
                jproc, 4, io_globalComm, status, ierr)
           
           call MPI_RECV(boundBoxt(1,1,1), 2*MDIM, & 
                FLASH_REAL, jproc, 6, io_globalComm, status, ierr)
           
           
           call MPI_RECV(gidt(1,1), ngid, & 
                FLASH_INTEGER, jproc, 7, io_globalComm, status, ierr)
           


!!******************************************************************************
!! OUTPUT MASTER_PE data.  No need to send this.
!!******************************************************************************
        else !end io_globalMe == MASTER_PE
      
           gidt(:,:)             = gr_gid(:,:)
           blockCenterCoordst(:,1) = blockCenterCoords(:)
           boundBoxt(:,:,1)        = boundBox(:,:)
           
        end if
   

        procnumber(1) = jproc
        call io_xfer_tree_data(tree_data, fileID, IO_FILE_HDF5, &
             IO_WRITE_XFER_MASTER_PE, &
             localNumBlocks, offset, presentDims)


        ! increment the offset 
        ! we just wrote localNumBlocks blocks from
        ! processor jproc to the output file, with the UG, localNumBlocks
        ! is always just 1
        offset = offset + localNumBlocks


     else ! if (io_globalMe == MasterPE)

        if (jproc == io_globalMe) then

                       
           call MPI_SEND(blockCenterCoords(1), MDIM, FLASH_REAL, &
                0, 4, io_globalComm, ierr)

           call MPI_SEND(boundBox(1,1), 2*MDIM, &
                FLASH_REAL, 0, 6, io_globalComm, ierr)

           call MPI_SEND(gr_gid(1,1), ngid, & 
                FLASH_INTEGER, 0, 7, io_globalComm, ierr)
           
        end if !! jproc == io_globalMe

     end if                 ! if io_globalMe == MASTER_PE



     !------------------------------------------------
     ! end processor loop
     !------------------------------------------------
  end do

!!Unk-type vars output faster if you loop over the variable first:  Moved them to here.

  ! Each unknown is stored in a separate record.  Loop over them, and pass the
  ! index into the he unk array, and the label of the current variable.  Note,
  ! a pointer to unkBuf in its entirety is passed -- this is already 
  ! contiguous.

!!******************************************************************************
!!UNK Variables: MASTERPE side of things:
!!******************************************************************************

  do i = UNK_VARS_BEGIN, UNK_VARS_END
     offset = 0 !offset must be reset
     allocate(unkt(NXB, NYB, NZB, 1))                               
     do jproc = 0, io_globalNumProcs-1

        unkt(:,:,:,:) = unk(i, io_ilo:io_ihi, io_jlo:io_jhi, io_klo:io_khi, :)

      
        !if(jproc /= MASTER_PE) then       
        if(io_globalMe == MASTER_PE .and. jproc /= MASTER_PE) then
           
           call MPI_RECV(unkt(1,1,1,1), &
                NXB*NYB*NZB, &
                FLASH_REAL, &
                jproc, 7+i, io_globalComm, &
                status, ierr)
           
        end if
        
        if (jproc == io_globalMe .and. io_globalMe /= MASTER_PE) then !send to MASTERPE
           
           call MPI_SEND( &
                unkt(1, 1, 1, 1), &
                NXB*NYB*NZB, &
                FLASH_REAL, &
                MASTER_PE, &
                7+i, &
                io_globalComm, &
                ierr)
           
        end if
                
        if(io_globalMe == MASTER_PE) then
           if(io_doublePrecision) then !writing to a checkpoint file
              call io_h5write_unknowns(io_globalMe, &
                   fileID, & 
                   NXB, & 
                   NYB, & 
                   NZB, & 
                   unkt(:,:,:,:), & 
                   globalVarMin(i), &
                   globalVarMax(i), &
                   io_unklabels(i), &
                   localNumBlocks, &
                   globalNumBlocks,  & 
                   offset)
              
           else
              
              call io_isPlotVar(i, isPlotVar,MAPBLOCK_UNK)
              
              if (isPlotVar) then
                 
                 if (io_plotfileGridQuantityDP) then
                    call io_h5write_unknowns(io_globalMe, &
                         fileID, & 
                         NXB, & 
                         NYB, & 
                         NZB, & 
                         unkt(:,:,:,:), & 
                         globalVarMin(i), &
                         globalVarMax(i), &
                         io_unklabels(i), &
                         localNumBlocks, &
                         globalNumBlocks,  & 
                         offset)
                 else
                    singleUnk(1:NXB,1:NYB,1:NZB,1) = & 
                         real(unkt(1:NXB,1:NYB,1:NZB,1), &
                         kind = single)
                    spMin = real(globalVarMin(i), kind = single)
                    spMax = real(globalVarMax(i), kind = single)
                    call io_h5write_unknowns_sp(io_globalMe, &
                         fileID, & 
                         NXB,   & 
                         NYB,   & 
                         NZB, & 
                         spMin, &
                         spMax, &
                         singleUnk(:,:,:,:),          & 
                         io_unklabels(i),  & 
                         localNumBlocks,  & 
                         globalNumBlocks,  & 
                         offset)
                 end if
                 
              endif !if plotVar
              
           end if !if io_doublePrecision
        end if !if io_globalMe == MASTER_PE
        
        offset = offset + 1 !we're in UG, only one block per proc.
        
     end do
     deallocate(unkt)

  end do !end unkvars
  
!!******************************************************************************
!!SCRATCH Variables: MASTERPE side of things:
!!******************************************************************************       
  
  allocate (scratcht(NXB, NYB, NZB, localNumBlocks))  
 
 
  do i= SCRATCH_GRID_VARS_BEGIN,SCRATCH_GRID_VARS_END
     !before we send, we should check to see if these are actually supposed to be plotted
     call io_isPlotVar(i, isPlotVar, MAPBLOCK_SCRATCH)
     if(.not.isPlotVar) CYCLE !move on to the next variable if we aren't outputting this scratch variable
  
     offset = 0 !reset offset.
     
     scratcht(:,:,:,:) = scratch(i,io_ilo:io_ihi, io_jlo:io_jhi, &
          io_klo:io_khi,1:localNumBlocks)
     
     do jproc = 0, io_globalNumProcs-1
        
        
        if(io_globalMe == MASTER_PE .and. jproc /= MASTER_PE) then
           
           call MPI_RECV(scratcht(1,1,1,1), &
                NXB*NYB*NZB*localNumBlocks, &
                FLASH_REAL, &
                jproc, 9+NUNK_VARS+i, io_globalComm, &
                status, ierr)
        end if
        if(io_globalMe == jproc .and. jproc /= MASTER_PE) then
           
           call MPI_SEND( &
                scratcht(1, 1, 1, 1), &
                NXB*NYB*NZB*localNumBlocks, &
                FLASH_REAL, &
                MASTER_PE, &
                9+NUNK_VARS+i, &
                io_globalComm, &
                ierr)
           
        end if
     
        
        !write the scratch grid vars if the user defines any in flash.par
        !we can use the same routine as when writing the unknowns.
        
        if (io_globalMe == MASTER_PE) then 
           call io_isPlotVar(i, isPlotVar, MAPBLOCK_SCRATCH)
           if(isPlotVar) then
              
              if(io_doublePrecision .or. io_plotfileGridQuantityDP) then
                 
                 call io_h5write_unknowns(io_globalMe, &
                      fileID, & 
                      NXB, & 
                      NYB, & 
                      NZB, & 
                      scratcht(:,:,:,1:localNumBlocks), & 
                      globalVarMinScratch(i), &
                      globalVarMaxScratch(i), &
                      io_plotGridVarStr(i), &
                      localNumBlocks, &
                      globalNumBlocks,  & 
                      offset)
                 
              else
                 
                 singleScratch(1:NXB,1:NYB,1:NZB,1:localNumBlocks) = & 
                      real(scratcht(:,:,:,1:localNumBlocks), kind = single)
                 
                 spMin = real(globalVarMinScratch(i), kind = single)                 
                 spMax = real(globalVarMaxScratch(i), kind = single)
                 
                 
                 call io_h5write_unknowns_sp(io_globalMe, &
                      fileID, & 
                      NXB,   & 
                      NYB,   & 
                      NZB, & 
                      spMin, &
                      spMax, &
                      singleScratch(:,:,:,1:localNumBlocks),          & 
                      io_plotGridVarStr(i),  & 
                      localNumBlocks,  & 
                      globalNumBlocks,  & 
                      offset)
                 
              end if !if io_doublePrecision
           end if !if isPlotVar
        end if !io_globalMe == MASTER_PE

        offset = offset + 1 !increment: we're in UG   
     end do
     deallocate(scratcht)
     
     
  end do  !end scratch
  
  !!******************************************************************************
  !!Face Centered Variables: MASTERPE side of things:
  !!******************************************************************************        
#if(NFACE_VARS > 0)
  
  allocate(facext(NXB+1,NYB,NZB,localNumBlocks))
  if(NDIM .gt. 1) allocate(faceyt(NXB,NYB+1,NZB,localNumBlocks))
  if(NDIM .gt. 2) allocate(facezt(NXB,NYB,NZB+1,localNumBlocks))
  
  do i = 1, NFACE_VARS

     offset = 0
     facext(:,:,:,:) = facevarx(i, io_ilo:io_ihi+1, io_jlo:io_jhi, io_klo:io_khi, :)
     if(NDIM .gt. 1) faceyt(:,:,:,:) = facevary(i, io_ilo:io_ihi, io_jlo:io_jhi+1, io_klo:io_khi, :)
     if(NDIM .gt. 2) facezt(:,:,:,:) = facevarz(i, io_ilo:io_ihi, io_jlo:io_jhi, io_klo:io_khi+1, :)
                                

     do jproc = 0, io_globalNumProcs-1
        
        if(io_globalMe == MASTER_PE .AND. jproc /= MASTER_PE) then
           
           call MPI_RECV(facext(1,1,1,1), &
                (NXB+1)*NYB*NZB, &
                FLASH_REAL, &
                jproc, 9+i+NUNK_VARS+NSCRATCH_GRID_VARS,&
                io_globalComm, &
                status, ierr)
           
           
           
           if(NDIM .gt. 1) then
              
              call MPI_RECV(faceyt(1,1,1,1), &
                   NXB*(NYB+1)*NZB, &
                   FLASH_REAL, &
                   jproc, 9+i+NUNK_VARS+NSCRATCH_GRID_VARS+NFACE_VARS, &
                   io_globalComm, &
                   status, ierr)
              
              
           end if !end NDIM .gt. 1
           if(NDIM .gt. 2) then
              
              call MPI_RECV(facezt(1,1,1,1), &
                   NXB*NYB*(NZB+1), &
                   FLASH_REAL, &
                   jproc, 9+i+NUNK_VARS+NSCRATCH_GRID_VARS+(NFACE_VARS*2),&
                   io_globalComm, &
                   status, ierr)
              
              
           end if  !end NDIM .gt. 2
        end if
        if (jproc == io_globalMe .and. jproc /= MASTER_PE) then
              
              
!!!POST FACEVAR SENDS
#if(NFACE_VARS > 0)
              
           call MPI_SEND( &
                facext(1, 1, 1, 1), &
                (NXB+1)*NYB*NZB, &
                FLASH_REAL, &
                MASTER_PE, &
                9+i+NUNK_VARS+NSCRATCH_GRID_VARS, &
                io_globalComm, &
                ierr)
           
              if(NDIM .gt. 1) then
                 call MPI_SEND( faceyt(1,1,1,1),&
                      NXB*(NYB+1)*NZB, &
                      FLASH_REAL, &
                      MASTER_PE, &
                      9+i+NUNK_VARS+NSCRATCH_GRID_VARS+NFACE_VARS, &
                      io_globalComm, &
                      ierr)
              end if
              if(NDIM .gt. 2) then
                 call MPI_SEND(facezt(1,1,1,1), &
                      NXB*NYB*(NZB+1), &
                      FLASH_REAL, &
                      MASTER_PE, &
                      9+i+NUNK_VARS+NSCRATCH_GRID_VARS+(NFACE_VARS*2), &
                      io_globalComm, &
                      ierr)
              end if
              
#endif
           end if
           
                           
        if(io_globalMe == MASTER_PE) then
           !write data
#if(NFACE_VARS > 0) 

           
           if(io_doublePrecision) then
              
              call io_h5write_unknowns(io_globalMe, &
                   fileID, &
                   NXB + 1, &
                   NYB, &
                   NZB, &
                   facext(:,:,:,:), &
                   globalFaceXMin(i), &
                   globalFaceXMax(i), &
                   io_faceXVarLabels(i), &
                   localNumBlocks, &
                   globalNumBlocks, &
                   offset)
              
              if(NDIM .gt. 1) then
                 
                 call io_h5write_unknowns(io_globalMe, &
                      fileID, &
                      NXB, &
                      NYB + 1, &
                      NZB, &
                      faceyt(:,:,:,:), &
                      globalFaceYMin(i), &
                      globalFaceYMax(i), &
                      io_faceYVarLabels(i), &
                      localNumBlocks, &
                      globalNumBlocks, &
                      offset)
                 
              end if
              if(NDIM .gt. 2) then
                 
                 call io_h5write_unknowns(io_globalMe, &
                      fileID, &
                      NXB, &
                      NYB, &
                      NZB + 1, &
                      facezt(:,:,:,:), &
                      globalFaceZMin(i), &
                      globalFaceZMax(i), &
                      io_faceZVarLabels(i), &
                      localNumBlocks, &
                      globalNumBlocks, &
                      offset)
                 
              end if
              
           end if !end io_doublePrecision
#endif
           
        end if !end io_globalMe == MASTER_PE
           
        offset = offset + 1   
     end do
  end do
  
  deallocate(facext)
  if(NDIM .gt. 1) deallocate(faceyt)
  if(NDIM .gt. 2) deallocate(facezt)
  
#endif  !! end NFACE_VARS > 0
  
  
  deallocate(gidt)
call MPI_BARRIER (io_globalComm, ierr)

  return
end subroutine io_writeData




