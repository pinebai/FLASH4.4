!!****if* source/Simulation/SimulationMain/unitTest/Multipole/IO_outputFinal
!!
!! NAME
!!
!!  IO_outputFinal
!!
!! SYNOPSIS
!!
!!  IO_outputFinal ()
!!
!! DESCRIPTION
!!
!!  This routine is called after the code has exited the main timestep loop.
!!  It outputs the last checkpoint, plotfile and particle plotfile.
!!  If particles are not included a stub (empty) routine will be called.
!!
!! ARGUMENTS
!!
!! NOTES
!!
!!  This is the overridden version for testing the grid multipole solver.  
!!  This ensures output to files.
!!
!!***

subroutine IO_outputFinal ()
  
  use IO_interface, ONLY : IO_writeCheckpoint, &
                           IO_writePlotfile,   &
                           IO_writeParticles

  implicit none

  call IO_writeCheckpoint ()
  call IO_writePlotfile   (.true.)
  call IO_writeParticles  (.false.)

end subroutine IO_outputFinal
