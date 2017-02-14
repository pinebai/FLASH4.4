!!****if* source/IO/IOMain/hdf5/parallel/PM/io_readData
!!
!! NAME
!!
!!  io_readData
!!
!!
!! SYNOPSIS
!!
!!  io_readData() 
!!
!!
!! DESCRIPTION
!!
!!  This is the reading counterpart to io_writeData.  It reads an HDF5
!!  file and distributes it to the processors to restart a simulation.
!!
!!  All reading is done using parallel HDF5 -- no explicit data movement
!!  is used.
!!
!!
!! ARGUMENTS
!!
!!***

!!REORDER(5): unk, facevar[xyz], unkBuf, unkBufGC, face[XYZ]Buf

#ifdef DEBUG_ALL
#define DEBUG_IO
#endif

#include "constants.h"
#include "Flash.h"
#include "io_flash.h"

subroutine io_readData()



  use Grid_data, ONLY : gr_globalNumBlocks, gr_nToLeft, &
       gr_globalOffset, gr_gid
  use Simulation_interface, ONLY : Simulation_mapStrToInt
  use Driver_interface, ONLY : Driver_abortFlash
  use Logfile_interface, ONLY : Logfile_stamp
  use Grid_interface, ONLY : Grid_putLocalNumBlks, Grid_getBlkIndexLimits, &
       Grid_receiveInputData

  use IO_data, ONLY : io_globalMe, io_globalComm,&
       io_baseName, io_checkpointFileNumber, &
       io_realParmNames, io_realParmValues, io_numRealParms, &
       io_intParmNames, io_intParmValues, io_numIntParms, &
       io_logParmNames, io_logParmValues, io_numLogParms, &
       io_strParmNames, io_strParmValues, io_numStrParms, &
       io_realScalarNames, io_realScalarValues, io_numRealScalars, &
       io_intScalarNames, io_intScalarValues, io_numIntScalars, &
       io_logScalarNames, io_logScalarValues, io_numLogScalars, &
       io_strScalarNames, io_strScalarValues, io_numStrScalars, &
       io_logToIntScalarValues, io_logToIntParmValues, io_unklabels, &
       io_ilo, io_ihi, io_jlo, io_jhi, io_klo, io_khi, &
       io_outputSplitNum, io_comm, io_chkptFileID, &
       io_chkGuardCellsInput, io_chkGuardCellsOutput, &
       io_faceXVarLabels, io_faceYVarLabels, io_faceZVarLabels,&
       io_realParmNamesPrev, io_realParmValuesPrev, io_numRealParmsPrev, &
       io_intParmNamesPrev, io_intParmValuesPrev, io_numIntParmsPrev, &
       io_logParmNamesPrev, io_logParmValuesPrev, io_numLogParmsPrev, &
       io_strParmNamesPrev, io_strParmValuesPrev, io_numStrParmsPrev, &
       io_logToIntParmValuesPrev, io_splitNumBlks, &
       io_unklabelsGlobal, io_meshMe, io_meshNumProcs, tree_data_t
  use IO_interface, ONLY : IO_getScalar

#ifdef FLASH_GRID_PARAMESH
  use tree, ONLY : maxblocks_tr, nfaces, nchild, nodetype, &
       lrefine, bnd_box, coord, bsize, neigh, parent, child
#ifdef FLASH_GRID_PARAMESH3OR4
  use tree, ONLY : MFLAGS, which_child, bflags
  use Grid_data, ONLY : gr_gsurr_blks
#endif
#endif
  use io_typeInterface, ONLY : io_xfer_tree_data
  use physicaldata, ONLY : unk, facevarx, facevary, facevarz


  implicit none
  type(tree_data_t) :: tree_data

#include "Flash_mpi.h"





  integer :: localNumBlocks, ngid
  integer :: alocalNumBlocks

 
  character (len=4) :: fnumStr
  character (len=MAX_STRING_LENGTH) :: filename
  character(len=MAX_STRING_LENGTH), save, allocatable, dimension(:,:) :: strBuff


  integer :: blockID, procBlocks, ierr

  integer :: i, lb, j, xx, yy, alnblocks, u

  integer, allocatable :: procnumber(:) !dimension(localNumBlocks)
  
  integer :: realGlobalNumBlocks
!  real :: unkBufGC(1,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC,MAXBLOCKS)
  real,allocatable :: unkBufGC(:,:,:,:,:)
  real,allocatable :: unkBuf(:,:,:,:,:)
  real,allocatable :: faceXBuf(:,:,:,:,:)
  real,allocatable :: faceYBuf(:,:,:,:,:)
  real,allocatable :: faceZBuf(:,:,:,:,:)
  
  integer :: blkLimits(2,MDIM), blkLimitsGC(2, MDIM)

  integer :: splitOffset, localOffset
  integer :: doread

  integer :: fileID, presentDims
  integer, parameter :: xferType = IO_READ_XFER, libType = IO_FILE_HDF5


  call io_getOutputName(io_checkpointFileNumber, "hdf5", "_chk_", filename, .false.)


  call io_h5open_file_for_read(io_chkptFileID, filename, io_comm, io_outputSplitNum)

  if (io_globalMe == MASTER_PE) then
       allocate (strBuff(2,2))
       print *, 'file: ', trim(filename), ' opened for restart'
       write (strBuff(1,1), "(A)") "type"
       write (strBuff(1,2), "(A)") "checkpoint"
       write (strBuff(2,1), "(A)") "name"
       write (strBuff(2,2), "(A)") trim(filename)
       call Logfile_stamp( strBuff, 2, 2, "[io_readData]")
  end if
  if (allocated(strBuff)) deallocate(strBuff)


 !!read in the total number of blocks, time, and timestep
  call io_h5read_header(io_globalMe, &
                        io_chkptFileID, &
                        io_unklabels, &
                        io_outputSplitNum)



  call io_prepareListsRead()
  


  call io_h5read_runtime_parameters(io_chkptFileID, &
       io_numRealParmsPrev, &
       io_realParmNamesPrev, &
       io_realParmValuesPrev, &
       io_numIntParmsPrev, &
       io_intParmNamesPrev, &
       io_intParmValuesPrev, &
       io_numStrParmsPrev, &
       io_strParmNamesPrev, &
       io_strParmValuesPrev, &
       io_numLogParmsPrev, &
       io_logParmNamesPrev, &
       io_logToIntParmValuesPrev)



  call io_h5read_scalars(io_chkptFileID, &
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
       io_logToIntScalarValues)


  call io_finalizeListsRead()


  call IO_getScalar("globalNumBlocks", gr_globalNumBlocks)
  call io_checkBlockShape(gr_globalNumBlocks)
  if (io_outputSplitNum /= 1) then
     call IO_getScalar("splitNumBlocks", io_splitNumBlks)
  else
     io_splitNumBlks = gr_globalNumBlocks
  endif


  !---------------------------------------------------------------------------
  ! compute the number of blocks on each processor -- this will be used to
  ! get the offset into the file for the parallel read
  !---------------------------------------------------------------------------

  ! compute the approximate number of blocks per processor
  alnblocks = int(gr_globalNumBlocks/io_meshNumProcs) + 1
  
  ! check for error -- if the number of blocks we want to put on each
  ! processor is greater than maxblocks, then abort
  if (alnblocks .GT. MAXBLOCKS) then
     
     print *
     print *, '********** ERROR in READ_DATA ************'
     print *
     print *,' Number of blocks per processor exceeds maxblocks.'
     print *,' Suggest you reset maxblocks to a larger number or'
     print *,' run on a larger number of processors.'
     print *,' globalNumBlocks, io_meshNumProcs = ', gr_globalNumBlocks, io_meshNumProcs
     print *
     
     call Driver_abortFlash('[io_readData] ERROR: num blocks per proc exceeds maxblocks')
     
  end if
  
  ! figure out the excess blocks
  yy = (io_meshNumProcs*alnblocks) - gr_globalNumBlocks
  xx = io_meshNumProcs - yy
  
  ! loop over all the processor numbers and figure out how many blocks are
  ! stored to the left of the processor -- this is a little tricky
  
  gr_nToLeft(0) = 0
  
  do i = 0, io_meshNumProcs - 2
     if (i .LT. xx) then
        procBlocks = alnblocks
     else
        procBlocks = alnblocks - 1
     endif
     
     if (alnblocks .EQ. 0) then
        if (i .LT. gr_globalNumBlocks) then
           procBlocks = 1
        else
           procBlocks = 0
        end if
     end if
     
     ! we have the number of blocks on proc i, the number of blocks on i+1 is
     ! the number of blocks on i + the number of blocks left of i
     if (i .EQ. 0) then
        gr_nToLeft(i+1) = procBlocks
     else
        gr_nToLeft(i+1) = procBlocks + gr_nToLeft(i)
     endif
  enddo
  
  ! figure out how many blocks are on the current proc.
  if (io_meshMe < xx) then
     localNumBlocks = alnblocks
  else
     localNumBlocks = alnblocks - 1
  endif
  
  if (alnblocks .EQ. 0) then
     if (io_meshMe < gr_globalNumBlocks) then
        localNumBlocks = 1
     else
        localNumBlocks = 0
     end if
  end if
  
  ! compute the offset into the dataspace in the HDF5 file
  gr_globalOffset = gr_nToLeft(io_meshMe)
  


  call Grid_putLocalNumBlks(localNumBlocks)

  !find our offset into a potentially split file:
  if(io_outputSplitNum > 1) then
     call MPI_ALLREDUCE(gr_globalOffset, splitOffset,1,FLASH_INTEGER,&
                        MPI_MIN, io_comm, ierr)
     localOffset = gr_globalOffset - splitOffset
     
  else
     localOffset = gr_globalOffset
  end if


  call io_h5read_present_dims(io_chkptFileID, presentDims)



  tree_data % bnd_box => bnd_box
  tree_data % coord => coord
  tree_data % bsize => bsize
  tree_data % gid => gr_gid
  tree_data % nodetype => nodetype
  tree_data % lrefine => lrefine
# ifdef FLASH_GRID_PARAMESH3OR4
  tree_data % bflags => bflags
  tree_data % which_child => which_child
# else
  nullify(tree_data % bflags)
  nullify(tree_data % which_child)
# endif
  nullify(tree_data % procnumber)
# ifdef FLASH_GRID_PARAMESH4DEV_SURR_BLKS_OPTION
  tree_data % gsurr_blks => gr_gsurr_blks
# else
  nullify(tree_data % gsurr_blks)
# endif

  call io_h5read_present_dims(io_chkptFileID, presentDims)

  call io_xfer_tree_data(tree_data, &
       io_chkptFileID, libType, xferType, &
       localNumBlocks, localOffset, presentDims)

  !Extract data from gr_gid and gr_gsurr_blks arrays
  call Grid_receiveInputData(localNumBlocks, alnblocks, xx)




  allocate(unkBuf(1, NXB, NYB, NZB, MAXBLOCKS))

  if (io_chkGuardCellsInput) then    ! currently, blkLimits is not used otherwise
     blockID = 1
     call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)
#ifndef FL_NON_PERMANENT_GUARDCELLS
     blkLimits = blkLimitsGC
#endif
  end if

  do u=1, ubound(io_unklabelsGlobal,1)
     call Simulation_mapStrToInt(io_unklabelsGlobal(u),i,MAPBLOCK_UNK)
     doread = 0
     if(i /= NONEXISTENT) doread = 1

     if(io_chkGuardCellsInput) then
        allocate(unkBufGC(1,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC,MAXBLOCKS))
        unkBufGC(1,:GRID_IHI_GC,:GRID_JHI_GC,:GRID_KHI_GC,1:localNumBlocks) = 0.0
        call io_h5read_unknowns(io_chkptFileID, &
             GRID_IHI_GC, &
             GRID_JHI_GC, &
             GRID_KHI_GC, &
             unkBufGC, &
             io_unklabelsGlobal(u), &
             localNumBlocks, &
             io_splitNumBlks,  &
             localOffset, &
             doread)
        
        if(i /= NONEXISTENT) then
           unk(i,:,:,:,1:MAXBLOCKS) = & 
                unkBufGC(1,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
                   blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS), &
                   blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS), 1:MAXBLOCKS)
        end if
        deallocate(unkBufGC)
     else
        unkBuf(1,:,:,:,1:localNumBlocks) = 0.0 
        call io_h5read_unknowns(io_chkptFileID, &
             NXB, &
             NYB, &
             NZB, &
             unkBuf, &
             io_unklabelsGlobal(u), &
             localNumBlocks, &
             io_splitNumBlks,  &
             localOffset, &
             doread)
        
        if(i /= NONEXISTENT) then
           unk(i,io_ilo:io_ihi, io_jlo:io_jhi, io_klo:io_khi,1:MAXBLOCKS) = & 
                unkBuf(1,1:NXB,1:NYB,1:NZB,1:MAXBLOCKS)
        end if        
     end if
  enddo

  deallocate(unkBuf)
#ifdef FLASH_GRID_PARAMESH3OR4
#if(NFACE_VARS>0)

    
    allocate(faceXBuf(1,1:NXB+1,1:NYB,1:NZB,1:MAXBLOCKS))
    if (NDIM .gt. 1) allocate(faceYBuf(1,1:NXB,1:NYB+1,1:NZB,1:MAXBLOCKS))
    if (NDIM .gt. 2) allocate(faceZBuf(1,1:NXB,1:NYB,1:NZB+1,1:MAXBLOCKS))

    doread = 1
    do i = 1,NFACE_VARS
      
      call io_h5read_unknowns(io_chkptFileID, &
                              NXB+1, &
                              NYB, &
                              NZB, &
                              faceXBuf, &
                              io_faceXVarLabels(i), &
                              localNumBlocks, &
                              io_splitNumBlks, &
                              localOffset, &
                              doread)
      
      facevarx(i,io_ilo:io_ihi+1,io_jlo:io_jhi,io_klo:io_khi,1:MAXBLOCKS) = &
         faceXBuf(1,1:NXB+1,1:NYB,1:NZB,1:MAXBLOCKS)
         
      if (NDIM .gt. 1) then
        
        call io_h5read_unknowns(io_chkptFileID, &
                                NXB, &
                                NYB+1, &
                                NZB, &
                                faceYBuf, &
                                io_faceYVarLabels(i), &
                                localNumBlocks, &
                                io_splitNumBlks, &
                                localOffset, &
                                doread)

        facevary(i,io_ilo:io_ihi,io_jlo:io_jhi+1,io_klo:io_khi,1:MAXBLOCKS) = &
          faceYBuf(1,1:NXB,1:NYB+1,1:NZB,1:MAXBLOCKS)
      
      end if

      if (NDIM .gt. 2) then
        call io_h5read_unknowns(io_chkptFileID, &
                                NXB, &
                                NYB, &
                                NZB+1, &
                                faceZBuf, &
                                io_faceZVarLabels(i), &
                                localNumBlocks, &
                                io_splitNumBlks, &
                                localOffset, &
                                doread)

        facevarz(i,io_ilo:io_ihi,io_jlo:io_jhi,io_klo:io_khi+1,1:MAXBLOCKS) = &
          faceZBuf(1,1:NXB,1:NYB,1:NZB+1,1:MAXBLOCKS)
       
       end if

    end do
    
    deallocate(faceXBuf)
    if (NDIM .gt. 1) deallocate(faceYBuf)
    if (NDIM .gt. 2) deallocate(faceZBuf)
    
#endif
#endif

  if (io_globalMe == MASTER_PE) &
    print *, 'read_data:  read ', gr_globalNumBlocks, ' blocks.'


  if (io_globalMe == MASTER_PE) then
     allocate (strBuff(2, 2))
     write (strBuff(1,1), "(A)") "type"
     write (strBuff(1,2), "(A)") "checkpoint"
     write (strBuff(2,1), "(A)") "name"
     write (strBuff(2,2), "(A)") trim(filename)
     call Logfile_stamp( strBuff, 2, 2, "[io_readData] file_closed")
     if (allocated(strBuff)) deallocate(strBuff)
  end if

 
  call MPI_BARRIER (io_globalComm, ierr)
  if (io_globalMe == MASTER_PE) &
    print *, 'io_readData:  finished reading input file.'

  return
end subroutine io_readData
