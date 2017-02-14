!!****if* source/IO/IOMain/hdf5/parallel/UG/io_readData
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
!! ARGUMENTS
!!
!!***

!!REORDER(5): unk, unkBuf, facevar[xyz], face[XYZ]Buf

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
       io_logToIntScalarValues, io_logToIntParmValuesPrev, &
       io_unklabels, io_unklabelsGlobal, &
       io_ilo, io_ihi, io_jlo, io_jhi, io_klo, io_khi, &
       io_checkpointFileNumber, io_baseName, io_comm, io_outputSplitNum, &
       io_chkptFileID, io_faceXVarLabels, io_faceYVarLabels, io_faceZVarLabels, &
       io_meshMe, io_meshNumProcs
  use Simulation_interface, ONLY : Simulation_mapStrToInt
  use Logfile_interface, ONLY : Logfile_stamp
  use Grid_interface, ONLY : Grid_putLocalNumBlks


  use Grid_data, ONLY : gr_globalOffset, gr_globalNumBlocks

  use physicalData, only : unk, facevarx, facevary, facevarz 

  implicit none

#include "Flash_mpi.h"
#include "constants.h"
#include "Flash.h"


  real, dimension(2, MDIM) :: boundBox, blkLimits, blkLimitsGC


  real :: blockCenterCoords(MDIM, 1) !!1 because only 1 block per proc
  real :: bsize(MDIM, 1) !!1 because only 1 block per proc

  integer :: localNumBlocks, ngid

  integer :: i, u, stat, ierr

  integer :: procnumber(1)

  integer :: globalNumBlocks, globalOffset

  
  character (len=4) :: fnumString
  character (len=MAX_STRING_LENGTH) :: filename
  character(len=MAX_STRING_LENGTH), save, allocatable, dimension(:,:) :: strBuff

  real,allocatable :: unkBuf(:,:,:,:,:)
  real,allocatable :: faceXBuf(:,:,:,:,:)
  real,allocatable :: faceYBuf(:,:,:,:,:)
  real,allocatable :: faceZBuf(:,:,:,:,:)

  integer :: doread


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

  call io_checkBlockShape(io_meshNumProcs)

  !if all procs are writing to the same file then globalNumBlocks = io_globalNumProcs
  gr_globalNumBlocks = io_meshNumProcs
  
  !if all procs are writing to the same file then globalOffset = io_globalMe
  gr_globalOffset = io_meshMe


  localNumBlocks = 1

  call Grid_putLocalNumBlks(localNumBlocks)


  allocate(unkBuf(1,NXB,NYB,NZB, 1))

  do u=1, ubound(io_unklabelsGlobal,1)
    call Simulation_mapStrToInt(io_unklabelsGlobal(u),i,MAPBLOCK_UNK)
    doread = 0
    if(i /= NONEXISTENT) doread = 1
    call io_h5read_unknowns(io_chkptFileID, &
          NXB, &
          NYB, &
          NZB, &
          unkBuf, &
          io_unklabelsGlobal(u), &
          localNumBlocks, &
          gr_globalNumBlocks,  &
          gr_globalOffset, &
          doread)

     if(i /= NONEXISTENT) then
       unk(i,io_ilo:io_ihi, io_jlo:io_jhi, io_klo:io_khi,1) = & 
            unkBuf(1,1:NXB,1:NYB,1:NZB, 1)
     end if
  enddo

  deallocate(unkBuf)

#if(NFACE_VARS>0)
    
    allocate(faceXBuf(1,NXB+1,NYB,NZB,1))
    if(NDIM .gt. 1) allocate(faceYBuf(1,1:NXB,1:NYB+1,1:NZB,1))
    if(NDIM .gt. 2) allocate(faceZBuf(1,1:NXB,1:NYB,1:NZB+1,1))
    
    do i = 1, NFACE_VARS
      doread = 1
      call io_h5read_unknowns(io_chkptFileID, &
                              NXB +1, &
                  NYB, &
                  NZB, &
                  faceXBuf, &
                  io_faceXVarLabels(i), &
                  localNumBlocks, &
                  gr_globalNumBlocks, &
                  gr_globalOffset, &
                  doread)
      
      facevarx(i,io_ilo:io_ihi+1,io_jlo:io_jhi,io_klo:io_khi,1) = &
        faceXBuf(1,1:NXB+1,1:NYB,1:NZB,1)
    
      if (NDIM .gt. 1) then

        call io_h5read_unknowns(io_chkptFileID, &
                              NXB, &
                  NYB + 1, &
                  NZB, &
                  faceYBuf, &
                  io_faceYVarLabels(i), &
                  localNumBlocks, &
                  gr_globalNumBlocks, &
                  gr_globalOffset, &
                  doread)
        facevary(i,io_ilo:io_ihi,io_jlo:io_jhi+1,io_klo:io_khi,1) = &
          faceYBuf(1,1:NXB,1:NYB+1,1:NZB,1)
      end if

      if (NDIM .gt. 2)  then

        call io_h5read_unknowns(io_chkptFileID, &
                              NXB, &
                  NYB, &
                  NZB + 1, &
                  faceZBuf, &
                  io_faceZVarLabels(i), &
                  localNumBlocks, &
                  gr_globalNumBlocks, &
                  gr_globalOffset, &
                  doread)
        facevarz(i,io_ilo:io_ihi,io_jlo:io_jhi,io_klo:io_khi+1,1) = &
          faceZBuf(1,1:NXB,1:NYB,1:NZB+1,1)
      end if
    end do

    deallocate(faceXBuf)
    if(NDIM .gt. 1) deallocate(faceYBuf)
    if(NDIM .gt. 2) deallocate(faceZBuf)
    
#endif


  call mpi_barrier (io_globalComm, ierr)
  if (io_globalMe == MASTER_PE) &
    print *, 'read_data:  read ', gr_globalNumBlocks, ' blocks.'



  return
end subroutine io_readData
