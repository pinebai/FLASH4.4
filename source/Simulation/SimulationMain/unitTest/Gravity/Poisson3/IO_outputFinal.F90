!!****if* source/Simulation/SimulationMain/unitTest/Gravity/Poisson3/IO_outputFinal
!!
!! NAME
!!
!!  IO_outputFinal
!!
!!
!! SYNOPSIS
!!
!!  IO_outputFinal() 
!!
!!
!! DESCRIPTION
!!
!!  This routine is called after the code has exited the main timestep
!!  loop.  It outputs the last checkpoint,plotfile and particle plotfile.
!!
!!  If particles are not included a stub (empty) routine will be called.
!!
!! ARGUMENTS
!!
!!
!!
!! NOTES
!!
!!  This is the overridden version for the Poisson3 unit test.  
!!  This overrides the io_justCheckpointed check to ensure that
!!  the final result is output.
!!
!!***

subroutine IO_outputFinal()
  
  use IO_data, ONLY : io_justCheckpointed
  use IO_interface, ONLY : IO_writeCheckpoint, IO_writePlotfile, &
    IO_writeParticles

  implicit none


  io_justCheckpointed = .false.
  if(.not. io_justCheckpointed) then
  
     call IO_writeCheckpoint()
     
     call IO_writePlotfile(.true.)
     
     call IO_writeParticles(.false.)
  end if

end subroutine IO_outputFinal
