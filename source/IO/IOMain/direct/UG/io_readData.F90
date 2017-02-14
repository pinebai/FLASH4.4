!!****if* source/IO/IOMain/direct/UG/io_readData
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
!!  This is the plain, non parallel, fortran read counterpart to io_writeData.  
!!  Each processor reads a binary file.
!!
!!  Order of reads matters!
!!
!!***


subroutine io_readData()
  
  use IO_data, ONLY : io_globalMe, io_globalNumProcs, io_realParmNamesPrev, io_realParmValuesPrev, &
       io_numRealParmsPrev, &
       io_intParmNamesPrev, io_intParmValuesPrev, io_numIntParmsPrev, &
       io_logParmNamesPrev, io_logParmValuesPrev, io_numLogParmsPrev, &
       io_strParmNamesPrev, io_strParmValuesPrev, io_numStrParmsPrev, &
       io_realScalarNames, io_realScalarValues, io_numRealScalars, &
       io_intScalarNames, io_intScalarValues, io_numIntScalars, &
       io_logScalarNames, io_logScalarValues, io_numLogScalars, &
       io_strScalarNames, io_strScalarValues, io_numStrScalars, &
       io_logToIntScalarValues, io_logToIntParmValuesPrev, io_unklabels, &
       io_ilo, io_ihi, io_jlo, io_jhi, io_klo, io_khi, &
       io_checkpointFileNumber, io_baseName, &
       io_chkptFileID
  use Logfile_interface, ONLY : Logfile_stamp
  use Grid_interface, ONLY : Grid_putLocalNumBlks


  use physicalData, only : unk 

  implicit none

#include "constants.h"
#include "Flash.h"

  real :: boundBox(2, MDIM)
  real :: blockCenterCoords(MDIM, 1) !!1 because only 1 block per proc
  real :: bsize(MDIM, 1) !!1 because only 1 block per proc


  integer :: localNumBlocks, ngid, ioLun, seedlength
  integer, allocatable :: fakegid(:,:)

  integer :: i, stat, ierr

  integer :: procnumber(1)

  integer :: fileID, alocalNumBlocks, globalNumBlocks, globalOffset

  integer :: blocksPerFile

  character (len=4) :: fnumString
  character (len=MAX_STRING_LENGTH) :: filename
  character(len=MAX_STRING_LENGTH), save, allocatable, dimension(:,:) :: strBuff


  call io_initFile(io_checkpointFileNumber, io_chkptFileID, filename, "_chk_", .false.)

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
  
  
  
  !first read the number of parameters and scalars
  read(io_chkptFileID)io_numRealParmsPrev
  read(io_chkptFileID)io_numIntParmsPrev
  read(io_chkptFileID)io_numStrParmsPrev
  read(io_chkptFileID)io_numLogParmsPrev
  
  read(io_chkptFileID)io_numRealScalars
  read(io_chkptFileID)io_numIntScalars
  read(io_chkptFileID)io_numStrScalars
  read(io_chkptFileID)io_numLogScalars


  !then allocate space for them
  call io_prepareListsRead()
  

  !! read the runtime parameters
  read(io_chkptFileID)io_realParmNamesPrev
  read(io_chkptFileID)io_realParmValuesPrev
  read(io_chkptFileID)io_intParmNamesPrev
  read(io_chkptFileID)io_intParmValuesPrev
  read(io_chkptFileID)io_strParmNamesPrev
  read(io_chkptFileID)io_strParmValuesPrev
  read(io_chkptFileID)io_logParmNamesPrev
  read(io_chkptFileID)io_logToIntParmValuesPrev
  
  !! read the scalars  
  read(io_chkptFileID)io_realScalarNames
  read(io_chkptFileID)io_realScalarValues
  read(io_chkptFileID)io_intScalarNames
  read(io_chkptFileID)io_intScalarValues
  read(io_chkptFileID)io_strScalarNames
  read(io_chkptFileID)io_strScalarValues
  read(io_chkptFileID)io_logScalarNames
  read(io_chkptFileID)io_logToIntScalarValues

  call io_finalizeListsRead()
  
  localNumBlocks = 1
  
  call Grid_putLocalNumBlks(localNumBlocks)

  if (NDIM == 1) then
     ngid = 5
  else if (NDIM == 2) then
     ngid = 9
  else if (NDIM == 3) then
     ngid = 15
  end if
  
  allocate(fakegid(ngid,1)) !1 because in UG only 1 block per proc      
  
  !these aren't used, but they are written out in checkpoint write
  !so we need to read them in
  read(io_chkptFileID)fakegid
  
  read(io_chkptFileID)boundBox
  
  read(io_chkptFileID)blockCenterCoords
  
  read(io_chkptFileID)bsize
  
  do i = UNK_VARS_BEGIN,UNK_VARS_END
     read(io_chkptFileID)unk(i,io_ilo:io_ihi, io_jlo:io_jhi, io_klo:io_khi,1) 
  enddo
  
  deallocate(fakegid)

  return

end subroutine io_readData
