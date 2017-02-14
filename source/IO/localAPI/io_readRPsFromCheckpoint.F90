!!****if* source/IO/localAPI/io_readRPsFromCheckpoint
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
!!  In the pnetcdf inplementation, all processors are involved in
!!  file operations even though only the information on MASTER_PE
!!  is ultimately required.
!!
!!
!!***
subroutine io_readRPsFromCheckpoint(myPE, filename, success, checkpointFileNumber)

  implicit none
  integer,intent(INOUT) :: myPE
  character(len=*),intent(IN) :: filename
  logical,intent(OUT) :: success
  integer,intent(OUT) :: checkpointFileNumber
  
  success = .FALSE.
  checkpointFileNumber = -1

end subroutine io_readRPsFromCheckpoint
