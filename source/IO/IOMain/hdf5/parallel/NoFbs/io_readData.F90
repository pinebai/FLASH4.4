!!****if* source/IO/IOMain/hdf5/parallel/NoFbs/io_readData
!!
!! NAME
!!
!!  io_readData
!!
!!
!! SYNOPSIS
!!
!!  call io_readData()
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
!! ARGUMENTS
!!   
!!   none
!!
!!***


!!REORDER(5): unk, facevar[xyz], unkBuf, face[XYZ]Buf

subroutine io_readData()
  
  use IO_data, ONLY : io_globalMe, io_globalNumProcs,io_globalComm,&
        io_realParmNamesPrev, io_realParmValuesPrev, io_numRealParmsPrev, &
       io_intParmNamesPrev, io_intParmValuesPrev, io_numIntParmsPrev, &
       io_logParmNamesPrev, io_logParmValuesPrev, io_numLogParmsPrev, &
       io_strParmNamesPrev, io_strParmValuesPrev, io_numStrParmsPrev, &
       io_realScalarNames, io_realScalarValues, io_numRealScalars, &
       io_intScalarNames, io_intScalarValues, io_numIntScalars, &
       io_logScalarNames, io_logScalarValues, io_numLogScalars, &
       io_strScalarNames, io_strScalarValues, io_numStrScalars, &
       io_logToIntScalarValues, io_logToIntParmValuesPrev, io_unklabels, io_unklabelsGlobal, &
       io_ilo, io_ihi, io_jlo, io_jhi, io_klo, io_khi, &
       io_checkpointFileNumber, io_baseName, &
       io_comm, io_outputSplitNum, io_chkptFileID,&
       io_iguard, io_jguard, io_kguard, &
       io_faceXVarLabels, io_faceYVarLabels, io_faceZVarLabels
  use Simulation_interface, ONLY : Simulation_mapStrToInt
  use Logfile_interface, ONLY : Logfile_stamp
  use Grid_interface, ONLY : Grid_putLocalNumBlks, &
    Grid_getBlkIndexLimits, Grid_getBlkCornerID
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get


  use physicalData, only : unk, facevarx, facevary, facevarz 

  implicit none

#include "Flash_mpi.h"
#include "constants.h"
#include "Flash.h"

  integer, dimension(LOW:HIGH, MDIM) :: boundBox
  integer, dimension(LOW:HIGH, MDIM) :: blkLimits, blkLimitsGC


  real :: blockCenterCoords(MDIM, 1) !!1 because only 1 block per proc
  real :: bsize(MDIM, 1) !!1 because only 1 block per proc

  integer :: localNumBlocks, ngid

  integer :: i, u, stat, ierr

  integer :: procnumber(1)

  integer :: globalNumBlocks, globalOffset

  integer :: blocksPerFile
  integer :: doread
  
  character (len=4) :: fnumString
  character (len=MAX_STRING_LENGTH) :: filename
  character(len=MAX_STRING_LENGTH), save, allocatable, dimension(:,:) :: strBuff

  real,allocatable :: unkBuf(:,:,:,:,:)
  real,allocatable :: faceXBuf(:,:,:,:,:)
  real,allocatable :: faceYBuf(:,:,:,:,:)
  real,allocatable :: faceZBuf(:,:,:,:,:)

  
  integer :: cornerID(MDIM), stride(MDIM)
  integer :: globalIndexLimits(MDIM)
  integer :: localnxb, localnyb, localnzb
  integer :: nxbOffset, nybOffset, nzbOffset, global_offset
  integer :: iprocs,jprocs,kprocs,iGridSize,jGridSize,kGridSize




  call io_getOutputName(io_checkpointFileNumber, "hdf5", "_chk_", filename, .false.)

  call io_h5open_file_for_read(io_chkptFileID, filename, io_comm, io_outputSplitNum)

  if (io_globalMe == MASTER_PE) then
       allocate (strBuff(2,2))
       print *, 'file: ', trim(filename), ' opened for restart'
       write (strBuff(1,1), "(A)") "type"
       write (strBuff(1,2), "(A)") "checkpoint"
       write (strBuff(2,1), "(A)") "name"
       write (strBuff(2,2), "(A)") trim(filename)
       call Logfile_stamp( strBuff, 2, 2, "[io_readData] file opened")
  end if
  if (allocated(strBuff)) deallocate(strBuff)



 !!read in header info
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


 
  !! now set all the parameters read in from the checkpoint to the "previous value"
  !! in the linked list dataset
  call io_finalizeListsRead()


  !find the number of procs (in UG, procs=blocks) writing to each split file
  !if (as is typical) only one file is written per dump meaning io_outputSplitNum = 1
  !procsPerFile simply  = io_globalNumProcs
  blocksPerFile = io_globalNumProcs / io_outputSplitNum

  !if all procs are writing to the same file then globalNumBlocks = io_globalNumProcs
  globalNumBlocks = blocksPerFile
  
  !call Grid_getGlobalIndexLimits(globalIndexLimits)
  
  !if all procs are writing to the same file then globalOffset = io_globalMe
  !but if splitting files then consider io_outputSplitNum
  globalOffset = mod(io_globalMe, blocksPerFile)
  globalOffset = 0

  localNumBlocks = 1

  call Grid_putLocalNumBlks(localNumBlocks)

  !This is non-fixed block size, NXB NYB and NZB mean nothing.
  !Need local versions instead.
  call Grid_getBlkIndexLimits(1, blkLimits,blkLimitsGC, CENTER)
  localnxb = blkLimits(HIGH,IAXIS) - blkLimits(LOW,IAXIS) + 1
  localnyb = blkLimits(HIGH,JAXIS) - blkLimits(LOW,JAXIS) + 1
  localnzb = blkLimits(HIGH,KAXIS) - blkLimits(LOW,KAXIS) + 1
  
  !This hack is to get around the fact that gr_cornerID is not set until 
  !Grid_initDomain is called.  We have to manually decompose the processors and 
  !figure out the offsets from a series of runtime parameters.

  call RuntimeParameters_get("iprocs", iprocs)
  call RuntimeParameters_get("jprocs", jprocs)
  call RuntimeParameters_get("kprocs", kprocs)
  call RuntimeParameters_get("iGridSize", iGridSize)
  call RuntimeParameters_get("jGridSize", jGridSize)
  call RuntimeParameters_get("kGridSize", kGridSize)

  !generate local nxb nyb and nzb
  localnxb = iGridSize / iprocs 
  localnyb = jGridSize / jprocs
  localnzb = kGridSize / kprocs
  
  !in the case that these don't divide evenly:
  if((mod(iGridsize, iprocs) .NE. 0) .AND. (io_globalMe .LT. mod(iGridSize, iprocs))) localnxb = localnxb + 1
  if((mod(jGridsize, jprocs) .NE. 0) .AND. (io_globalMe .LT. mod(jGridSize, jprocs))) localnxb = localnxb + 1
  if((mod(kGridsize, kprocs) .NE. 0) .AND. (io_globalMe .LT. mod(kGridSize, kprocs))) localnxb = localnxb + 1
    

  call Grid_getBlkCornerID(1,cornerID,stride)
  
  nxbOffset = cornerID(IAXIS) - 1
  nybOffset = (cornerID(JAXIS) - 1) * K2D  
  nzbOffset = (cornerID(KAXIS) - 1) * K3D

!!$  nxbOffset = (mod(io_globalMe, iprocs) *localnxb)
!!$  nybOffset = (mod((io_globalMe / iprocs), jprocs) * localnyb)
!!$  nzbOffset = (mod((io_globalMe / (iprocs * jprocs)), kprocs) * localnzb)
!!$  
!!$  !if (nxbOffset .ne. 0) nxbOffset = nxbOffset + 1
!!$  !if (nybOffset .ne. 0) nybOffset = nybOffset + 1
!!$  !if (nzbOffset .ne. 0) nzbOffset = nzbOffset  1

  globalIndexLimits(IAXIS) = iGridSize
  globalIndexLimits(JAXIS) = jGridSize
  globalIndexLimits(KAXIS) = kGridSize
  
  print *, "localnb values are ", localnxb,localnyb,localnzb
  print *, "offsets are ", nxbOffset,nybOffset,nzbOffset
  print *, "io_guard", io_iguard, io_jguard, io_kguard

  allocate(unkBuf(1,localnxb,localnyb,localnzb, 1))
  
  !do i = UNK_VARS_BEGIN,UNK_VARS_END
  do u=1, ubound(io_unklabelsGlobal,1)
    call Simulation_mapStrToInt(io_unklabelsGlobal(u),i,MAPBLOCK_UNK)
    doread = 0
    if(i /= NONEXISTENT) doread = 1
    call io_h5read_unknowns(io_chkptFileID, &
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
      io_unklabelsGlobal(u), &
      doread)
    
    if(i /= NONEXISTENT) then
      unk(i,1+io_iguard:localnxb+io_iguard,&
            1+io_jguard:localnyb+io_jguard,&
            1+io_kguard:localnzb+io_kguard,1) = &
            unkBuf(1,1:localnxb,1:localnyb,1:localnzb,1)
    end if
  enddo

  deallocate(unkBuf)

#if(NFACE_VARS > 0)

  allocate(faceXBuf(1,localnxb+1,localnyb,localnzb,1))
  if(NDIM > 1) allocate(faceYBuf(1,localnxb,localnyb+1,localnzb,1))
  if(NDIM > 2) allocate(faceZBuf(1,localnxb,localnyb,localnzb+1,1))

  do i = 1, NFACE_VARS

    call io_h5read_facevars(io_chkptFileID, &
           globalIndexLimits(IAXIS) + 1, &
           globalIndexLimits(JAXIS), &
           globalIndexLimits(KAXIS), &
           nxbOffset, &
           nybOffset, &
           nzbOffset, &
           localnxb +1, &
           localnyb, &
           localnzb, &
           faceXBuf, &
           io_faceXVarLabels(i))
           
     facevarx(i, 1+io_iguard:localnxb+io_iguard+1, &
                 1+io_jguard:localnyb+io_jguard, &
                 1+io_kguard:localnzb+io_kguard, 1) = &
       faceXBuf(1,1:localnxb+1,1:localnyb,1:localnzb,1)

    if( NDIM > 1) then
     
      call io_h5read_facevars(io_chkptFileID, &
             globalIndexLimits(IAXIS), &
             globalIndexLimits(JAXIS) + 1, &
             globalIndexLimits(KAXIS), &
             nxbOffset, &
             nybOffset, &
             nzbOffset, &
             localnxb, &
             localnyb + 1, &
             localnzb, &
             faceYBuf, &
             io_faceYVarLabels(i))
           
       facevary(i, 1+io_iguard:localnxb+io_iguard, &
                   1+io_jguard:localnyb+io_jguard+1, &
                   1+io_kguard:localnzb+io_kguard, 1) = &
         faceYBuf(1,1:localnxb,1:localnyb+1,1:localnzb,1)

    end if !NDIM > 1


    if(NDIM > 2) then
    
      call io_h5read_facevars(io_chkptFileID, &
             globalIndexLimits(IAXIS), &
             globalIndexLimits(JAXIS), &
             globalIndexLimits(KAXIS) + 1, &
             nxbOffset, &
             nybOffset, &
             nzbOffset, &
             localnxb, &
             localnyb, &
             localnzb + 1, &
             faceZBuf, &
             io_faceZVarLabels(i))
           
       facevarz(i, 1+io_iguard:localnxb+io_iguard, &
                   1+io_jguard:localnyb+io_jguard, &
                   1+io_kguard:localnzb+io_kguard+1, 1) = &
         faceZBuf(1,1:localnxb,1:localnyb,1:localnzb+1,1)

    end if !NDIM > 2

  end do

  deallocate(faceXBuf)
  if(NDIM > 1) deallocate(faceYBuf)
  if(NDIM > 2) deallocate(faceZBuf)
 
#endif


  call mpi_barrier (io_globalComm, ierr)
  if (io_globalMe == MASTER_PE) &
    print *, 'read_data:  read ', globalNumBlocks, ' blocks.'



  return
end subroutine io_readData
