!!****if* source/IO/IOMain/IO_sendOutputData
!!
!! NAME
!!  IO_sendOutputData
!!
!! SYNOPSIS
!! 
!!  IO_sendOutputData()
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
!!  IO_updateScalars
!!
!!***

subroutine IO_sendOutputData()

  use IO_data, ONLY : io_doublePrecision, io_nextCheckpointTime, io_nextPlotfileTime, io_splitNumBlks, io_splitParts, &
       io_nextPlotFileZ, io_nextCheckpointZ, &
       io_plotFileNumber, io_forcedPlotFileNumber, io_checkpointFileNumber
  use IO_interface, ONLY : IO_setScalar
  implicit none

!!$  print *, "inside IO_sendOutputData"
  
  !CD: "corners" does nothing.  However, it is part of the checkpoint file 
  !specification, so we retain the checkpoint entry, but force it to be .false..
  !It was used in FLASH2 to interpolate the data to the zone corners before 
  !storing the data in the plotfile (for creating improved iso-surfaces).
  call IO_setScalar("corners", .false.)
  call IO_setScalar("double_precision", io_doublePrecision)
  call IO_setScalar("nextCheckpointTime", io_nextCheckpointTime)
  call IO_setScalar("nextPlotfileTime", io_nextPlotfileTime)
  call IO_setScalar("nextCheckpointZ", io_nextCheckpointZ)
  call IO_setScalar("nextPlotFileZ", io_nextPlotFileZ)
  
  !used for split IO
  call IO_setScalar("splitNumBlocks", io_splitNumBlks)
  
  call IO_setScalar("splitNumParticles", io_splitParts)
  
  !used for restarting with -c checkpointFile without a flash.par
  call IO_setScalar('checkpointFileNumber', io_checkpointFileNumber)
  call IO_setScalar('plotFileNumber', io_plotFileNumber)
  call IO_setScalar('forcedPlotFileNumber', io_forcedPlotFileNumber)

  !we may bave particles included
  call io_ptSendOutputData()

end subroutine IO_sendOutputData

