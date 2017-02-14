!!****if* source/IO/IOMain/hdf5/serial/PM/io_writeData
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

!!REORDER(5): unk, scratch, facevar[xyz], scratchBuf

#include "constants.h"
#include "Flash.h"
#include "io_flash.h"

subroutine io_writeData (fileID)

  use io_typeInterface, ONLY : io_xfer_tree_data

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
       io_ilo, io_ihi, io_jlo, io_jhi, io_klo, io_khi, io_geometry, &
       io_setupCall, io_buildDir, io_flashRelease, &
       io_fileCreationTime, io_buildDate, io_buildMachine, io_cflags, io_fflags, &
       io_setupTimeStamp, io_buildTimeStamp, io_outputSplitNum, io_doublePrecision, &
       io_plotVarStr, io_nPlotVars, io_plotGridVarStr, io_scratchGridVarlabels,&
       io_faceXVarLabels, io_faceYVarLabels, io_faceZVarLabels, &
       io_plotFaceVarStr, io_plotfileMetadataDP, io_plotfileGridQuantityDP, &
       io_fileFormatVersion, tree_data_t
  use Grid_interface, ONLY : Grid_getLocalNumBlks

  use Grid_data, ONLY : gr_globalNumBlocks, gr_gid, scratch

#ifdef FLASH_GRID_PARAMESH
  use tree, ONLY : maxblocks_tr, nfaces, nchild, nodetype, &
       lrefine, bnd_box, coord, bsize
#ifdef FLASH_GRID_PARAMESH3OR4
  use Grid_data, ONLY : gr_gsurr_blks
  use tree, ONLY : MFLAGS, which_child, bflags
#endif
#endif

  use physicaldata, ONLY :unk,facevarx,facevary,facevarz

  implicit none


#include "Flash_mpi.h"


  integer, INTENT(in) :: fileID

  type(tree_data_t) :: tree_data
  integer :: jproc, i, j, blockID
  integer :: localNumBlockst 
  integer :: ierr


  ! Define the temp variables (t) which hold the data transfered
  ! to the MASTER_PE
  integer, target :: lrefinet(MAXBLOCKS), nodetypet(MAXBLOCKS)
  integer, target :: gidt(nfaces+1+nchild,MAXBLOCKS)
  integer :: status(MPI_STATUS_SIZE)

  integer,allocatable :: procnumber(:)

  real, target :: coordt(MDIM, MAXBLOCKS), sizet(MDIM,MAXBLOCKS)
  real, target :: bnd_boxt(2,MDIM, MAXBLOCKS) 
  real, allocatable :: unkt(:,:,:,:)
  real, allocatable :: scratcht(:,:,:,:)
  real, allocatable :: faceXt(:,:,:,:)
  real, allocatable :: faceYt(:,:,:,:)
  real, allocatable :: faceZt(:,:,:,:)

  real :: scratchBuf(NSCRATCH_GRID_VARS, NXB, NYB, NZB, MAXBLOCKS)


  integer :: localNumBlocks

  ! storage for the global block number we are writing out
  integer :: offset

  ! for logfile output
  character(len=16)                                    :: numToStr


  ! block data message buffering stuff. not using this right now
  integer, parameter :: MAX_TRANS_SIZE = 4000000  ! Maximum number of bytes

  ! to transfer at one time
  ! when sending block data
  integer, save      :: MAX_BLK_MSGS
  integer            :: n_blk_msgs, pblkcount


  ! allocate storage to hold a single variable information
  ! this should only be a small memory overhead
  integer, parameter :: single = SELECTED_REAL_KIND(p=6)

  !  real (kind=single) :: singleUnk(1,NXB,NYB,NZB,MAXBLOCKS)
  real (kind=single), allocatable :: singleUnk(:,:,:,:)
  real (kind=single) :: singleScratch(NXB,NYB,NZB,MAXBLOCKS)
  !Face Centered:
  real (kind=single) :: singleFaceX(NXB+1,NYB,NZB,MAXBLOCKS)
  real (kind=single) :: singleFaceY(NXB,NYB+1,NZB,MAXBLOCKS)
  real (kind=single) :: singleFaceZ(NXB,NYB,NZB+1,MAXBLOCKS)


  ! allocate storage to hold the coordinate information and bounding box
  ! information
  real (kind=single) :: sizeSingle(MDIM,MAXBLOCKS)
  real (kind=single) :: coordSingle(MDIM,MAXBLOCKS)
  real (kind=single) :: bndSingle(2,MDIM,MAXBLOCKS)
  real (kind=single) :: spMax, spMin

  logical :: isPlotVar
  logical, allocatable :: isScratchPlotVar(:)
  real  :: globalVarMin(NUNK_VARS), globalVarMax(NUNK_VARS)
  real  :: globalVarMinScratch(NSCRATCH_GRID_VARS), globalVarMaxScratch(NSCRATCH_GRID_VARS)
  real :: globalVarMinFaceX(NFACE_VARS), globalVarMaxFaceX(NFACE_VARS)
  real :: globalVarMinFaceY(NFACE_VARS), globalVarMaxFaceY(NFACE_VARS)
  real :: globalVarMinFaceZ(NFACE_VARS), globalVarMaxFaceZ(NFACE_VARS)


#ifdef FLASH_GRID_PARAMESH3OR4
  !data structures specific to PM3 f
  integer, target :: bflagst(mflags, MAXBLOCKS), which_childt(MAXBLOCKS), &
       gsurr_blkst(2,1+(K1D*2),1+(K2D*2),1+(K3D*2),MAXBLOCKS)
#endif
  integer, parameter :: presentDims = MDIM, libType = IO_FILE_HDF5

  if (io_globalMe == MASTER_PE) then
     call io_createDatasets(fileID, gr_globalNumBlocks, presentDims)
  end if

  tree_data % bnd_box => bnd_boxt
  tree_data % coord => coordt
  tree_data % bsize => sizet
  tree_data % gid => gidt
  tree_data % nodetype => nodetypet
  tree_data % lrefine => lrefinet
#ifdef FLASH_GRID_PARAMESH3OR4
  tree_data % bflags => bflagst
  tree_data % which_child => which_childt
  tree_data % gsurr_blks => gsurr_blkst
#else
  nullify(tree_data % bflags)
  nullify(tree_data % which_child)
  nullify(tree_data % gsurr_blks)
#endif
  allocate(tree_data % procnumber(MAXBLOCKS))


  call Grid_getLocalNumBlks(localNumBlocks)

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

     ! keep track of the block written out in a number of processors independent
     ! manner
     offset = 0

  end if


  !get the max and minimum variables
  call io_getVarExtrema(NUNK_VARS, globalVarMin, globalVarMax, CENTER)

  allocate(isScratchPlotVar(SCRATCH_GRID_VARS_BEGIN:SCRATCH_GRID_VARS_END))
  do i = SCRATCH_GRID_VARS_BEGIN,SCRATCH_GRID_VARS_END
     call io_isPlotVar(i, isScratchPlotVar(i), MAPBLOCK_SCRATCH)
  end do
  if (ANY(isScratchPlotVar)) then
     call io_getVarExtrema(NSCRATCH_GRID_VARS, globalVarMinScratch, globalVarMaxScratch, SCRATCH)
  end if

#ifdef FLASH_GRID_PARAMESH3OR4
#if(NFACE_VARS > 0) 
  call io_getVarExtrema(NFACE_VARS, globalVarMinFaceX, globalVarMaxFaceX, FACEX)

  if(NDIM .gt. 1) &
       call io_getVarExtrema(NFACE_VARS, globalVarMinFaceY, globalVarMaxFaceY, FACEY)

  if(NDIM .gt. 2) &
       call io_getVarExtrema(NFACE_VARS, globalVarMinFaceZ, globalVarMaxFaceZ, FACEZ)
#endif  
#endif

  !-----------------------------------------------------------------------------
  ! loop over all of the processors.  All the data is moved to processor 0 for
  ! storage using MPI sends and receives.
  !-----------------------------------------------------------------------------
  localNumBlockst=0
  do jproc = 0,io_globalNumProcs-1

     if (io_globalMe == MASTER_PE) then

        ! fetch localNumblocks from other processors
        if (jproc /= 0) then

           call MPI_RECV (localNumBlockst,1,FLASH_INTEGER,jproc, & 
                1,io_globalComm,status,ierr)

           if (localNumBlockst > 0) then

              call MPI_RECV(lrefinet(1), localNumBlockst, FLASH_INTEGER, & 
                   jproc, 2, io_globalComm, status, ierr)

              call MPI_RECV(nodetypet(1),localNumBlockst, FLASH_INTEGER, & 
                   jproc, 3, io_globalComm, status, ierr)

              call MPI_RECV(coordt(1,1), MDIM*localNumBlockst, FLASH_REAL, &
                   jproc, 4, io_globalComm, status, ierr)

              call MPI_RECV(sizet(1,1), MDIM*localNumBlockst, FLASH_REAL, &
                   jproc, 5, io_globalComm, status, ierr)

              call MPI_RECV(bnd_boxt(1,1,1), 2*MDIM*localNumBlockst, & 
                   FLASH_REAL, jproc, 6, io_globalComm, status, ierr)

              call MPI_RECV(gidt(1,1), localNumBlockst*(nfaces+1+nchild), & 
                   FLASH_INTEGER, jproc, 7, io_globalComm, status, ierr)

#ifdef FLASH_GRID_PARAMESH3OR4
              call MPI_RECV(which_childt(1), localNumBlockst, FLASH_INTEGER, & 
                   jproc, 8, io_globalComm, status, ierr)

              call MPI_RECV(bflagst(1,1), localNumBlockst*MFLAGS, & 
                   FLASH_INTEGER, jproc, 9, io_globalComm, status, ierr)

              call MPI_RECV(gsurr_blkst(1,1,1,1,1), &
                   2*(1+(K1D*2))*(1+(K2D*2))*(1+(K3D*2))*localNumBlockst, &
                   FLASH_INTEGER, jproc, 9, io_globalComm, status, ierr)
#endif


              if(ANY(isScratchPlotVar)) then
                 allocate (scratcht(NXB, NYB, NZB, localNumBlockst))                               
                 do i= SCRATCH_GRID_VARS_BEGIN,SCRATCH_GRID_VARS_END

                    call MPI_RECV(scratcht(1,1,1,1), &
                         NXB*NYB*NZB*localNumBlockst, &
                         FLASH_REAL, &
                         jproc, 9+NUNK_VARS+i, io_globalComm, &
                         status, ierr)

                    scratchBuf(i,:,:,:,1:localNumBlockst) = scratcht(:,:,:,:)

                 end do
                 deallocate(scratcht)
              end if


           end if !(if localNumBlockst > 0)

        else

           localNumBlockst = localNumBlocks
           if (localNumBlockst > 0) then
              lrefinet(1:localNumBlockst)     = lrefine(1:localNumBlockst)
              nodetypet(1:localNumBlockst)    = nodetype(1:localNumBlockst)
              gidt(:,1:localNumBlockst)       = gr_gid(:,1:localNumBlockst)

              bnd_boxt(:,:,1:localNumBlockst) = bnd_box(:,:,1:localNumBlockst)
              coordt(:,1:localNumBlockst)     = coord(:,1:localNumBlockst)
              sizet(:,1:localNumBlockst)      = bsize(:,1:localNumBlockst)

#ifdef FLASH_GRID_PARAMESH3OR4
              !paramesh3 specific data structures
              which_childt(1:localNumBlockst) = which_child(1:localNumBlockst)
              bflagst(:,1:localNumBlockst)    = bflags(:,1:localNumBlockst)
              gsurr_blkst(:,:,:,:,1:localNumBlockst) = &
                   gr_gsurr_blks(:,:,:,:,1:localNumBlockst)
#endif

              if(ANY(isScratchPlotVar)) then
                 scratchBuf(SCRATCH_GRID_VARS_BEGIN:SCRATCH_GRID_VARS_END,&
                      1:NXB,1:NYB,1:NZB,&
                      1:localNumBlockst)= & 
                      scratch(SCRATCH_GRID_VARS_BEGIN:SCRATCH_GRID_VARS_END,&
                      io_ilo:io_ihi,io_jlo:io_jhi,io_klo:io_khi,&
                      1:localNumBlockst)
              end if
           end if
        end if


        ! Write tree and mesh data
        tree_data % procnumber(:) = jproc
        call io_xfer_tree_data(tree_data, fileID, IO_FILE_HDF5, &
             IO_WRITE_XFER_MASTER_PE, &
             localNumBlockst, offset, presentDims)


        if (localNumBlockst > 0) then
           ! Each unknown is stored in a separate record.  Loop over them, and pass the
           ! index into the unk array, and the label of the current variable.  Note,
           ! a pointer to unkBuf in its entirety is passed -- this is already 
           ! contiguous.

           !deallocate(unkBuf)

           !write the scratch grid vars if the user defines any in flash.par
           !we can use the same routine as when writing the unknowns.
           if(ANY(isScratchPlotVar)) then
              do i = SCRATCH_GRID_VARS_BEGIN,SCRATCH_GRID_VARS_END

                 isPlotVar = isScratchPlotVar(i)
                 if(isPlotVar) then
                    if(io_doublePrecision .or. io_plotfileGridQuantityDP) then

                       call io_h5write_unknowns(io_globalMe, &
                            fileID, & 
                            NXB, & 
                            NYB, & 
                            NZB, & 
                            scratchBuf(i,:,:,:,1:localNumBlockst), & 
                            globalVarMinScratch(i), &
                            globalVarMaxScratch(i), &
                            io_scratchGridVarlabels(i), &
                            localNumBlockst, &
                            gr_globalNumBlocks,  & 
                            offset)

                    else

                       singleScratch(1:NXB,1:NYB,1:NZB,1:localNumBlockst) = & 
                            real(scratchBuf(i,:,:,:,1:localNumBlockst), kind = single)

                       spMin = real(globalVarMinScratch(i), kind = single)
                       spMax = real(globalVarMaxScratch(i), kind = single)



                       call io_h5write_unknowns_sp(io_globalMe, &
                            fileID, & 
                            NXB,   & 
                            NYB,   & 
                            NZB, & 
                            spMin, &
                            spMax, &
                            singleScratch(:,:,:,1:localNumBlockst),          & 
                            io_scratchGridVarlabels(i),  & 
                            localNumBlockst,  & 
                            gr_globalNumBlocks,  & 
                            offset)

                    end if !if io_doublePrecision
                 end if !if isPlotVar
              end do  !SCRATCH vars loop
           end if                !if ANY(isScratchPlotVar)

#ifdef FLASH_GRID_PARAMESH3OR4
           !Paramesh2 is incompatable with face-centered variables.
           !start outputing facevars if we have them.

#endif



        end if !if localNumBlockst > 0



        !--------------------------------------------------------------------
        ! end local block loop
        !--------------------------------------------------------------------

        ! increment the global block number -- 
        !we just wrote localNumBlockst blocks from
        ! processor jproc to the output file
        offset = offset + localNumBlockst


     else ! if (io_globalMe == MasterPE)

        if (jproc == io_globalMe) then

           call MPI_SEND(localNumblocks, 1, FLASH_INTEGER, 0, & 
                1, io_globalComm, ierr)

           if (localNumBlocks > 0) then

              call MPI_SEND(lrefine(1), localNumBlocks, FLASH_INTEGER, 0, & 
                   2, io_globalComm, ierr)

              call MPI_SEND(nodetype(1), localNumBlocks, FLASH_INTEGER, 0, & 
                   3, io_globalComm, ierr)

              call MPI_SEND(coord(1,1), MDIM*localNumBlocks, FLASH_REAL, &
                   0, 4, io_globalComm, ierr)

              call MPI_SEND(bsize(1,1), MDIM*localNumBlocks, FLASH_REAL, &
                   0, 5, io_globalComm,ierr)

              call MPI_SEND(bnd_box(1,1,1), 2*MDIM*localNumBlocks, &
                   FLASH_REAL, 0, 6, io_globalComm, ierr)

              call MPI_SEND(gr_gid(1,1), localNumBlocks*(nfaces+1+nchild), & 
                   FLASH_INTEGER, 0, 7, io_globalComm, ierr)

#ifdef FLASH_GRID_PARAMESH3OR4
              call MPI_SEND(which_child(1), localNumBlocks, FLASH_INTEGER, 0, & 
                   8, io_globalComm, ierr)

              call MPI_SEND(bflags(1,1), localNumBlocks*MFLAGS, & 
                   FLASH_INTEGER, 0, 9, io_globalComm, ierr)

              call MPI_SEND(gr_gsurr_blks(1,1,1,1,1), &
                   2*(1+(K1D*2))*(1+(K2D*2))*(1+(K3D*2))*localNumBlocks, &
                   FLASH_INTEGER, 0, 9, io_globalComm, ierr)
#endif



              if(ANY(isScratchPlotVar)) then
                 allocate(scratcht(NXB, NYB, NZB, localNumBlocks))
                 do i= SCRATCH_GRID_VARS_BEGIN,SCRATCH_GRID_VARS_END

                    scratcht(:,:,:,:) = scratch(i, io_ilo:io_ihi, io_jlo:io_jhi, &
                         io_klo:io_khi, 1:localNumBlocks)

                    call MPI_SEND( &
                         scratcht(1, 1, 1, 1), &
                         NXB*NYB*NZB*localNumBlocks, &
                         FLASH_REAL, &
                         MASTER_PE, &
                         9+NUNK_VARS+i, &
                         io_globalComm, &
                         ierr)
                 end do
                 deallocate(scratcht)
              end if

           endif !if localNumBlocks > 0

        end if !(if jproc == io_globalMe

     end if ! if io_globalMe == MASTER_PE



     !------------------------------------------------------------------------
     ! end processor loop
     !------------------------------------------------------------------------
  end do
  deallocate(isScratchPlotVar)

  !!*****************************************************************************
  !!output unk variables
  !!*****************************************************************************

  do i= UNK_VARS_BEGIN, UNK_VARS_END

     offset = 0
     do jproc = 0, io_globalNumProcs-1

        !allocate(unkBuf(NUNK_VARS, NXB, NYB, NZB, localNumBlockst))




        if(io_globalMe == MASTER_PE ) then

           if(jproc /= MASTER_PE) then        

              !post recieves
              call MPI_RECV (localNumBlockst,1,FLASH_INTEGER,jproc, & 
                   1,io_globalComm,status,ierr)

              if (localNumBlockst > 0) then
                 allocate(unkt(NXB, NYB, NZB, localNumBlockst))                                          

                 call MPI_RECV(unkt, &
                      NXB*NYB*NZB*localNumBlockst, &
                      FLASH_REAL, &
                      jproc, 9+i, io_globalComm, &
                      status, ierr)
                 !write out immediately
              end if


           else !we are on MASTER_PE

              localNumBlockst = localNumBlocks

              if (localNumBlockst > 0) then
                 allocate(unkt(NXB,NYB,NZB, localNumBlockst))
                 unkt(:,:,:,:) = unk(i, io_ilo:io_ihi, io_jlo:io_jhi, &
                      io_klo:io_khi, 1:localNumBlocks)
              end if
           end if !end if jprocs /= MASTER_PE


           if (localNumBlockst > 0) then

              if(io_doublePrecision) then
                 call io_h5write_unknowns(io_globalMe, &
                      fileID, & 
                      NXB, & 
                      NYB, & 
                      NZB, & 
                      unkt(:,:,:,1:localNumBlockst), & 
                      globalVarMin(i), &
                      globalVarMax(i), &
                      io_unklabels(i), &
                      localNumBlockst, &
                      gr_globalNumBlocks,  & 
                      offset)

              else

                 call io_isPlotVar(i, isPlotVar, MAPBLOCK_UNK)

                 if (isPlotVar) then

                    if(io_plotfileGridQuantityDP) then
                       call io_h5write_unknowns(io_globalMe, &
                            fileID, & 
                            NXB, & 
                            NYB, & 
                            NZB, & 
                            unkt(:,:,:,1:localNumBlockst), & 
                            globalVarMin(i), &
                            globalVarMax(i), &
                            io_unklabels(i), &
                            localNumBlockst, &
                            gr_globalNumBlocks,  & 
                            offset)
                    else
                       allocate(singleUnk(NXB,NYB,NZB,localNumBlockst))
                       singleUnk(1:NXB,1:NYB,1:NZB,1:localNumBlockst) = &
                            real(unkt(:,:,:,1:localNumBlockst), kind = single)
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
                            localNumBlockst,  & 
                            gr_globalNumBlocks,  & 
                            offset)
                       deallocate(singleUnk)
                    end if

                 endif !if plotVar

                 !unkBuf(i,:,:,:,1:localNumBlockst) = unkt(:,:,:,:)
              end if !if io_doublePrecision

              offset = offset + localNumBlockst
              deallocate(unkt)
           end if !if localNumBlockst > 0
        end if !if io_globalMe == MASTER_PE

        if (jproc == io_globalMe .and. jproc /= MASTER_PE) then

           call MPI_SEND(localNumblocks, 1, FLASH_INTEGER, 0, & 
                1, io_globalComm, ierr)

           if (localNumBlocks > 0) then
              allocate(unkt(NXB, NYB, NZB, localNumBlocks))

              unkt(:,:,:,:) = unk(i, io_ilo:io_ihi, io_jlo:io_jhi, &
                   io_klo:io_khi, 1:localNumBlocks)

              call MPI_SEND( &
                   unkt, &
                   NXB*NYB*NZB*localNumBlocks, &
                   FLASH_REAL, &
                   MASTER_PE, &
                   9+i, &
                   io_globalComm, &
                   ierr)

              deallocate(unkt)
           end if


        end if !end jproc /=MASTER_PE ^ jproc == io_globalMe

     end do !end jproc

  end do !end nvars loop

  !!*****************************************************************************
  !!output scratch vars
  !!*****************************************************************************

  do i = SCRATCH_GRID_VARS_BEGIN, SCRATCH_GRID_VARS_END

  end do

  !!*****************************************************************************
  !!output facevars only if we have them
  !!*****************************************************************************

#if NFACE_VARS > 0


  do i=1,NFACE_VARS

     offset = 0

     do jproc = 0,io_globalNumProcs - 1



        if(MASTER_PE == io_globalMe .and. MASTER_PE /= jproc) then
           !POST RECVS
           call MPI_RECV (localNumBlockst,1,FLASH_INTEGER,jproc, & 
                1,io_globalComm,status,ierr)

           if (localNumBlockst > 0) then
              allocate(facext(NXB+1,NYB,NZB,localNumBlockst))
              if(NDIM .gt. 1) allocate(faceyt(NXB,NYB+1,NZB,localNumBlockst))
              if(NDIM .gt. 2) allocate(facezt(NXB,NYB,NZB+1,localNumBlockst))

              call MPI_RECV(facext(1,1,1,1),&
                   (NXB+1)*NYB*NZB*localNumBlockst, &
                   FLASH_REAL, &
                   jproc, &
                   9+NUNK_VARS+NSCRATCH_GRID_VARS+i,&
                   io_globalComm, status, ierr)


              if(NDIM .GT. 1) then

                 call MPI_RECV(faceyt(1,1,1,1),&
                      NXB*(NYB+1)*NZB*localNumBlockst, &
                      FLASH_REAL, &
                      jproc, &
                      9+NUNK_VARS+NSCRATCH_GRID_VARS+NFACE_VARS+i,&
                      io_globalComm, status, ierr)
              end if !NDIM .GT. 1

              if(NDIM .GT. 2) then

                 call MPI_RECV(facezt(1,1,1,1), &
                      NXB*NYB*(NZB+1)*localNumBlockst, &
                      FLASH_REAL, &
                      jproc, &
                      9+NUNK_VARS+NSCRATCH_GRID_VARS+(NFACE_VARS*2)+i,&
                      io_globalComm, status, ierr)

              end if !NDIM .GT. 2
           end if
        end if !end recvs if

        if(MASTER_PE /= io_globalMe .AND. jproc == io_globalMe) then

           !POST SENDS

           !First we send the localNumBlocks
           call MPI_SEND(localNumblocks, 1, FLASH_INTEGER, 0, & 
                1, io_globalComm, ierr)

           if (localNumBlocks > 0) then
              !pack facevars for sending:
              allocate(facext(NXB+1,NYB,NZB,localNumBlocks))
              if(NDIM .gt. 1) allocate(faceyt(NXB,NYB+1,NZB,localNumBlocks))
              if(NDIM .gt. 2) allocate(facezt(NXB,NYB,NZB+1,localNumBlocks))

              facext(:,:,:,:) = facevarx(i, io_ilo:io_ihi+1, io_jlo:io_jhi, io_klo:io_khi, 1:localNumBlocks)
              if(NDIM .gt. 1) faceyt(:,:,:,:) = facevary(i, io_ilo:io_ihi, io_jlo:io_jhi+1, io_klo:io_khi, 1:localNumBlocks)
              if(NDIM .gt. 2) facezt(:,:,:,:) = facevarz(i, io_ilo:io_ihi, io_jlo:io_jhi, io_klo:io_khi+1, 1:localNumBlocks)


              !send what facevars we have now
              call MPI_SEND( facext(1,1,1,1), &
                   (NXB+1)*NYB*NZB*localNumBlocks, &
                   FLASH_REAL, &
                   MASTER_PE, &
                   9+NUNK_VARS+NSCRATCH_GRID_VARS+i, &
                   io_globalComm, &
                   ierr)

              if(NDIM .GT. 1) then

                 call MPI_SEND( faceyt(1,1,1,1), &
                      NXB*(NYB+1)*NZB*localNumBlocks, &
                      FLASH_REAL, &
                      MASTER_PE, &
                      9+NUNK_VARS+NSCRATCH_GRID_VARS+NFACE_VARS+i, &
                      io_globalComm, &
                      ierr)
              end if

              if(NDIM .GT. 2) then

                 call MPI_SEND(facezt(1,1,1,1), &
                      NXB*NYB*(NZB+1)*localNumBlocks, &
                      FLASH_REAL, &
                      MASTER_PE, &
                      9+NUNK_VARS+NSCRATCH_GRID_VARS+(NFACE_VARS*2)+i, &
                      io_globalComm, &
                      ierr)
              end if !end 3d

              !clean up jprocs' memory
              deallocate(facext)
              if(NDIM .gt. 1) deallocate(faceyt)
              if(NDIM .gt. 2) deallocate(facezt)
           end if !if localNumBlocks > 0

        end if !end sends if

        if(MASTER_PE == io_globalMe) then

           if(io_globalMe == jproc) then
              !we have to pack our own face*t array

              if (localNumBlocks > 0) then
                 allocate(facext(NXB+1,NYB,NZB,localNumBlocks))
                 if(NDIM .gt. 1) allocate(faceyt(NXB,NYB+1,NZB,localNumBlocks))
                 if(NDIM .gt. 2) allocate(facezt(NXB,NYB,NZB+1,localNumBlocks))

                 facext(:,:,:,:) = facevarx(i, io_ilo:io_ihi+1, io_jlo:io_jhi, io_klo:io_khi, 1:localNumBlocks)
                 if(NDIM .gt. 1) faceyt(:,:,:,:) = facevary(i, io_ilo:io_ihi, io_jlo:io_jhi+1, io_klo:io_khi, 1:localNumBlocks)
                 if(NDIM .gt. 2) facezt(:,:,:,:) = facevarz(i, io_ilo:io_ihi, io_jlo:io_jhi, io_klo:io_khi+1, 1:localNumBlocks)
              end if

              localNumBlockst = localNumBlocks

           end if !end local packing

           if (localNumBlockst > 0) then
              !we can write now
              if(io_doublePrecision) then
                 !we must be at least 1d
                 call io_h5write_unknowns(io_globalMe, &
                      fileID, &
                      NXB + 1, &
                      NYB, & 
                      NZB, &
                      facext(:,:,:,1:localNumBlockst), &
                      globalVarMinFaceX(i), &
                      globalVarMaxFaceX(i), &
                      io_faceXVarLabels(i), &
                      localNumBlockst, &
                      gr_globalNumBlocks, &
                      offset)
                 !are we 2d?

                 if(NDIM .GT. 1) then

                    call io_h5write_unknowns(io_globalMe, &
                         fileID, &
                         NXB, &
                         NYB + 1, & 
                         NZB, &
                         faceyt(:,:,:,1:localNumBlockst), &
                         globalVarMinFaceY(i), &
                         globalVarMaxFaceY(i), &
                         io_faceYVarLabels(i), &
                         localNumBlockst, &
                         gr_globalNumBlocks, &
                         offset)
                    !are we 3d?
                    if(NDIM .GT. 2) then

                       call io_h5write_unknowns(io_globalMe, &
                            fileID, &
                            NXB, &
                            NYB, & 
                            NZB + 1, &
                            facezt(:,:,:,1:localNumBlockst), &
                            globalVarMinFaceZ(i), &
                            globalVarMaxFaceZ(i), &
                            io_faceZVarLabels(i), &
                            localNumBlockst, &
                            gr_globalNumBlocks, &
                            offset)
                    end if !3d
                 end if !2d


              else !we're outputing a plot file
                 !Plotfile output not implemented at this point for facevars.

              end if !end for if(io_doublePrecision)

              deallocate(facext)
              if(NDIM .gt. 1) deallocate(faceyt)
              if(NDIM .gt. 2) deallocate(facezt)
           end if
        end if !end write if

        offset = offset + localNumBlockst

     end do ! end do jproc

  end do ! end do NFACE_VARS

#endif


  deallocate(tree_data % procnumber)
  nullify(tree_data % procnumber)

  call MPI_BARRIER (io_globalComm, ierr)

  return
end subroutine io_writeData
