!!****if* source/IO/IOMain/hdf5/io_readRPsFromCheckpoint
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
!!  In the hdf5 inplementation, the MASTER_PE alone does all
!!  the work.
!!
!! HISTORY
!!  Created - Klaus Weide
!!***

subroutine io_readRPsFromCheckpoint(myPE, filename, success, checkpointFileNumber)

  use RuntimeParameters_interface, ONLY : RuntimeParameters_set
  use IO_data, ONLY : io_baseName, io_checkpointFileNumber, &
       io_realParmNames, io_realParmValues, io_numRealParms, &
       io_intParmNames, io_intParmValues, io_numIntParms, &
       io_logParmNames, io_logParmValues, io_numLogParms, &
       io_strParmNames, io_strParmValues, io_numStrParms, &
       io_logToIntScalarValues, io_logToIntParmValues, io_unklabels, &
       io_comm, io_chkptFileID, &
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
#include "Flash_mpi.h"
  integer,intent(INOUT) :: myPE
  character(len=*),intent(IN) :: filename
  logical,intent(OUT) :: success
  integer,intent(OUT) :: checkpointFileNumber
  character(len=MAX_STRING_LENGTH)          :: strValue
  integer :: ipos, spos, pos, tempNum, ioStatus, mySplitNum
  integer :: ierr

  if (io_comm == 0) io_comm = io_globalComm
  if (io_comm == 0) then
     io_comm = MPI_COMM_WORLD
     call MPI_Comm_rank(io_comm, myPE, ierr)
  end if

  if (myPE == MASTER_PE) then
     success = .TRUE.

     ipos = 1
     spos = 1
     do while (ipos > 0)
        ipos = index(filename(spos:), '/')
        if (ipos > 0) spos = spos + ipos
     end do
     pos = index(filename(spos:), 'hdf5_chk_')
     if (pos > 0) then
        pos = pos + 9
     else
        pos = len_trim(filename(spos:)) - 3
     end if
     pos = max(1,pos)
     strValue = filename(spos+pos-1:spos+pos+3)
     !!  read (strValue, fmt='(I4)', iostat=ioStatus) tempNum
     read (strValue, fmt=*, iostat=ioStatus) tempNum
     if (ioStatus == 0) then
        !!     io_checkpointFileNumber = tempNum
        checkpointFileNumber = tempNum
     else
        success = .FALSE.
        print*,'Could not extract a checkpointFileNumber from "',trim(filename),'"!'
        return
     end if

     !prepare the io communicator if necessary - IO_init probably has not been called yet.

     !!  call io_getOutputName(io_checkpointFileNumber, "hdf5", "_chk_", filename, .false.)

     mySplitNum = 1
     call io_h5open_file_for_read(io_chkptFileID, filename, io_comm, mySplitNum)


     call io_h5read_runtime_parameters(io_chkptFileID, &
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
          io_logToIntParmValues)

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

     call io_h5close_file(io_chkptFileID)

     !!  call RuntimeParameters_set('checkpointFileNumber', tempNum)
  else
     success = .TRUE.
     checkpointFileNumber = -1
  end if
  return
end subroutine io_readRPsFromCheckpoint
