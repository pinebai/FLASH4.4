!!****if* source/IO/localAPI/io_ptSendOutputData
!!
!! NAME
!!  io_ptSendOutputData
!!
!! SYNOPSIS
!! 
!!  io_ptSendOutputData()
!!  
!! DESCRIPTION
!!
!!  This routine adds current values of variables belonging to the IO unit that should be
!!  checkpointed to the list of variables that are written out.
!!
!!  This routine is normally called by IO_updateScalars.
!!
!! NOTES
!!
!!  Variables that begin with "io_" like io_cellDataType and
!!  io_precision are stored in the data Fortran module for the 
!!  IO unit, IO_data.  The "io_" is meant to
!!  indicate that the variable belongs to the IO Unit.
!!
!! SEE ALSO
!!
!!
!!***

subroutine io_ptSendOutputData()
  
  implicit none
  
  return
 
end subroutine io_ptSendOutputData
