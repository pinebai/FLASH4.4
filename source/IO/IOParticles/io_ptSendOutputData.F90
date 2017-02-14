!!****if* source/IO/IOParticles/io_ptSendOutputData
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
  
  use IOParticles_data, ONLY : io_nextParticleFileZ, io_nextParticleFileTime, io_particleFileNumber

  use IO_interface, ONLY : IO_setScalar

  implicit none
  
  call IO_setScalar('nextParticleFileZ', io_nextParticleFileZ)
  call IO_setScalar('nextParticleFileTime', io_nextParticleFileTime)

  !used for restarting with -c checkpointFile without a flash.par
  call IO_setScalar('particleFileNumber', io_particleFileNumber)

  return
 
end subroutine io_ptSendOutputData
