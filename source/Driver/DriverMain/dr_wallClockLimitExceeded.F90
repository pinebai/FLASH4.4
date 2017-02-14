!!****if* source/Driver/DriverMain/dr_wallClockLimitExceeded
!!
!! NAME
!!
!!  dr_wallClockLimitExceeded
!!
!! SYNOPSIS
!!
!!  call dr_wallClockLimitExceeded(logical(OUT)  :: endRunWallClock)
!!
!! DESCRIPTION
!!
!!  Determine whether the wall clock time limit has been exceeded.
!!
!! ARGUMENTS
!!
!!   endRunWallClock : the answer is returned here.
!!
!! NOTES
!!
!!  This routine must be called collectively by all MPI tasks.
!!  The elapsed time is determiend on the master PE and then
!!  broadcast to all tasks.
!!
!!  This routines is usually called from Driver_evolveFlash.
!!
!!***

subroutine dr_wallClockLimitExceeded(endRunWallClock)
  use Driver_interface, ONLY: Driver_getElapsedWCTime
  use Driver_data, ONLY: dr_globalMe, dr_globalComm, &
                         dr_elapsedWCTime, dr_wallClockTimeLimit
  implicit none

  logical, intent(OUT) :: endRunWallClock

#include "constants.h"
  include "Flash_mpi.h"
  integer :: ierr

  endRunWallClock = .false.
  if (dr_wallClockTimeLimit .GE. 0.0) then
     call Driver_getElapsedWCTime(dr_elapsedWCTime)
     if (dr_globalMe == MASTER_PE) then
        if (dr_elapsedWCTime >  dr_wallClockTimeLimit) endRunWallClock = .true.
     end if
     call MPI_Bcast(endRunWallClock, 1, FLASH_LOGICAL, MASTER_PE, dr_globalComm, ierr)
  end if

end subroutine dr_wallClockLimitExceeded
