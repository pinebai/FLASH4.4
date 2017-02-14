!!****if* source/Simulation/SimulationMain/unitTest/IO/Driver_evolveFlash
!!
!! NAME
!!
!!  Driver_evolveFlash
!!
!! SYNOPSIS
!!
!!  Driver_evolveFlash()
!!
!! DESCRIPTION
!!
!! This is a very simple version of the Driver_evolveFlash routine,
!! that is meant to be used exclusively with IO Unit
!! testing. There is no time advancement involved here.
!!
!!  
!!
!! NOTES
!!
!! The Driver unit uses a few unit scope variables that are
!! accessible to all routines within the unit, but not to the
!! routines outside the unit. These variables begin with "dr_"
!! like, dr_globalMe or dr_dt, dr_beginStep, and are stored in fortran
!! module Driver_data (in file Driver_data.F90. The other variables
!! are local to the specific routine and do not have the prefix "dr_"
!!
!!
!!***

subroutine Driver_evolveFlash()

  use Driver_data, ONLY:   dr_nbegin,  dr_restart, dr_initialSimTime
  use Logfile_interface, ONLY : Logfile_stamp, Logfile_close
  use Timers_interface, ONLY : Timers_start, Timers_stop, &
    Timers_getSummary

  use physicaldata, ONLY : unk
  use IO_interface, ONLY : IO_writeCheckpoint, IO_writePlotfile, &
                           IO_writeParticles
  implicit none

#include "constants.h"
#include "Flash.h"

  integer :: iOut, var

  iOut = 2

  open(iOut,file='unitTest_0000')

  !initialize the fake grid variable with dummy values
  do var = 1, NUNK_VARS
     
     unk(var, :,:,:,:) = var*1.0

  end do
     
  call Timers_start("IO_writeCheckpoint")
  call IO_writeCheckpoint()
  call Timers_stop("IO_writeCheckpoint")

  call Timers_start("IO_writePlotfile")
  call IO_writePlotfile()
  call Timers_stop("IO_writePlotfile")

  call Timers_start("IO_writeParticles")
  call IO_writeParticles( .false.)
  call Timers_stop("IO_writeParticles")


  call Timers_getSummary(0)


  call Logfile_stamp( "FLASH run complete.", "LOGFILE_END")

  call Logfile_close()


  write(iOut,'("all results conformed with expected values.")')
  close(iOut)

  return

end subroutine Driver_evolveFlash

