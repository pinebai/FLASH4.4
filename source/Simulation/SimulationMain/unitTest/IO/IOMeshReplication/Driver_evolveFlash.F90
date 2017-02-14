!!****if* source/Simulation/SimulationMain/unitTest/IO/IOMeshReplication/Driver_evolveFlash
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

  if (.not.dr_restart) then
     !initialize the fake grid variable with dummy values
     do var = 1, NUNK_VARS
        unk(var,:,:,:,:) = var*1.0
     end do
     call init_nonrep_vars()
  end if

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


subroutine init_nonrep_vars
  use physicaldata, ONLY : unk
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_data, ONLY : dr_meshAcrossMe, dr_meshCopyCount, dr_meshMe, &
        dr_meshNumProcs
  implicit none
  integer :: totalSharedVars, sharedVarsInMyMesh, n, &
       globalUnkIndex, localUnkIndex

  call RuntimeParameters_get("totalSharedVars", totalSharedVars)
  sharedVarsInMyMesh = &
       NONREP_NLOCS(dr_meshAcrossMe, dr_meshCopyCount, totalSharedVars)

  !We have already done a generic initialization of the variables in unk.
  !Now we will re-initialize the shared unk variables with a global index
  !rather than a local index.  This will allow us to obtain the same data
  !in the checkpoint file irrespective of the number of meshes.
  do n = 1, sharedVarsInMyMesh
     globalUnkIndex = NONREP_LOC2GLOB(n, dr_meshAcrossMe, dr_meshCopyCount)
     localUnkIndex = X_NONREP_LOC2UNK(n)
     if (dr_meshMe == MASTER_PE) then
        write(6,'(4(a,i5))') " Mesh master:", dr_meshAcrossMe, &
             " - my team of size", dr_meshNumProcs, &
             " will write to unk index", localUnkIndex, &
             " which is non-replicated variable", globalUnkIndex
     end if
     unk(localUnkIndex,:,:,:,:) = globalUnkIndex*1.0
  end do
end subroutine init_nonrep_vars
