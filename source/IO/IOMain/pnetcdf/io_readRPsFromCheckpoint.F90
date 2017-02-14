!!****if* source/IO/IOMain/pnetcdf/io_readRPsFromCheckpoint
!!
!! NAME
!!
!!  io_readRPsFromCheckpoint
!!
!! SYNOPSIS
!!
!!  call io_readRPsFromCheckpoint(integer(in)  :: mype,
!!                                character(len=*)(IN)  :: filename,
!!                                logical(OUT)  :: success,
!!                                integer(OUT)  :: checkpointfilenumber)
!!
!! DESCRIPTION
!!
!!   Read lists of runtime parameters and scalars from
!!   an HDF5 file written previously by FLASH.
!!
!!   On successful return, the master PE will have names and values
!!   of runtime parameters and scalars in certain arrays that
!!   belong to the IO_data module:
!!      io_{real,int,log,str}Parm{Names,Values} and
!!      io_{real,int,log,str}Scalar{Names,Values}.
!!   The numbers of slots availabe in this arrays (which must be
!!   allocated before calling this subroutine) are given on
!!   entry as
!!      io_num{Real,Int,Log,Str}Parms and
!!      io_num{Real,Int,Log,Str}Scalars.
!!   On successful return, these variables contain the numbers
!!   of valid data items read from the file.
!!
!!   Only the master PE need do anything if that is how the IO
!!   implementation works.
!!
!! ARGUMENTS
!!
!!   mype : my Processor Number
!!
!!   filename : names a file from which runtime parameter and scalar
!!              information should be read.
!!
!!   success : indicates whether the implementation successfully
!!             parsed the given filename for a checkpoint file number
!!             and then opened the file given by the filename and
!!             populated the io_* lists for parameter and scalar
!!             names and values.
!!
!!   checkpointfilenumber : the number extracted from the filename is
!!             returned here. Note that this does not know anything about
!!             the checkpoint file number that may be stored as a scalar
!!             IN the file.
!!
!! NOTES
!!
!!  In the pnetcdf inplementation, all processors are involved in
!!  file operations even though only the information on MASTER_PE
!!  is ultimately required.
!!
!!
!!***

subroutine io_readRPsFromCheckpoint(myPE, filename, success, checkpointFileNumber)

  use RuntimeParameters_interface, ONLY : RuntimeParameters_set
  use IO_data, ONLY : io_baseName, io_checkpointFileNumber, &
       io_maxParms, &
       io_realParmNames, io_realParmValues, io_numRealParms, &
       io_intParmNames, io_intParmValues, io_numIntParms, &
       io_logParmNames, io_logParmValues, io_numLogParms, &
       io_strParmNames, io_strParmValues, io_numStrParms, &
       io_logToIntScalarValues, io_logToIntParmValues, io_unklabels, &
       io_outputSplitNum, io_comm, io_chkptFileID, &
       io_faceXVarLabels, io_faceYVarLabels, io_faceZVarLabels,&
       io_realParmNamesPrev, io_realParmValuesPrev, io_numRealParmsPrev, &
       io_intParmNamesPrev, io_intParmValuesPrev, io_numIntParmsPrev, &
       io_logParmNamesPrev, io_logParmValuesPrev, io_numLogParmsPrev, &
       io_strParmNamesPrev, io_strParmValuesPrev, io_numStrParmsPrev, &
       io_realScalarNames, io_realScalarValues, io_numRealScalars, &
       io_intScalarNames, io_intScalarValues, io_numIntScalars, &
       io_logScalarNames, io_logScalarValues, io_numLogScalars, &
       io_strScalarNames, io_strScalarValues, io_numStrScalars, &
       io_logToIntScalarValues, &
       io_logToIntParmValuesPrev, io_globalComm
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
  integer,intent(IN) :: myPE
  character(len=*),intent(IN) :: filename
  logical,intent(OUT) :: success
  integer,intent(OUT) :: checkpointFileNumber
  character(len=MAX_STRING_LENGTH)          :: strValue, fnameBuffer
  integer :: ipos, spos, pos, tempNum, ioStatus, strtype, istat
  

  success = .TRUE.

  call MPI_Type_Contiguous (MAX_STRING_LENGTH, MPI_CHARACTER, strtype, istat)
  call MPI_Type_Commit (strtype, istat)
  if (myPE == MASTER_PE) then
     fnameBuffer = filename
  end if
  call MPI_Bcast (fnameBuffer, 1, strtype, MASTER_PE, io_globalComm, istat)
  call MPI_Type_Free (strtype, istat)
   

  ipos = 1
  spos = 1
  do while (ipos > 0)
     ipos = index(fnameBuffer(spos:), '/')
     if (ipos > 0) spos = spos + ipos
  end do
  pos = index(fnameBuffer(spos:), 'ncmpi_chk_')
  if (pos > 0) then
     pos = pos + 10
  else
     pos = len_trim(fnameBuffer(spos:)) - 3
  end if
  pos = max(1,pos)
  strValue = fnameBuffer(spos+pos-1:spos+pos+3)
!!  read (strValue, fmt='(I4)', iostat=ioStatus) tempNum
  read (strValue, fmt=*, iostat=ioStatus) tempNum
  if (ioStatus == 0) then
!!     io_checkpointFileNumber = tempNum
    checkpointFileNumber = tempNum
  else
     success = .FALSE.
     print*,'Could not extract a checkpointFileNumber from "',trim(fnameBuffer),'"!'
     return
  end if
#ifdef DEBUG_ALL
  print*,myPE,': got fnameBuffer=="',fnameBuffer,'" is',checkpointFileNumber
#endif

!!  call io_getOutputName(io_checkpointFileNumber, "ncmpi", "_chk_", fnameBuffer, .false.)

#ifdef DEBUG_ALL
  print*,myPE,': will io_ncmpi_open_file_for_read,fnameBuffer=="',fnameBuffer,'"'
#endif
  call mpi_barrier(io_globalComm,ioStatus)

  call io_ncmpi_open_file_for_read(io_chkptFileID, fnameBuffer)
#ifdef DEBUG_ALL
  print*,myPE,': did io_ncmpi_open_file_for_read,fnameBuffer=="',fnameBuffer,'"', io_chkptFileID
#endif

  call io_ncmpi_read_header(myPE, &
       io_chkptFileID, &
       NUNK_VARS, &
       io_unklabels, &
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

  call io_ncmpi_close_file(io_chkptFileID)

!!  call RuntimeParameters_set('checkpointFileNumber', tempNum)

  return
end subroutine io_readRPsFromCheckpoint
