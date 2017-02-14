!!****if* source/IO/IOMain/io_restrictBeforeWrite
!!
!! NAME
!!  io_restrictBeforeWrite
!!
!! SYNOPSIS
!!
!!  io_restrictBeforeWrite() logical(OUT) :: restrictNeeded)
!!
!!
!! DESCRIPTION
!!
!!  Determines whether grid data must be restricted up the 
!!  tree before we perform a write.  
!!
!! ARGUMENTS
!!
!!  restrictNeeded : whether a restrict is needed.
!!
!! NOTES
!! 
!!  A restriction ensures all levels have valid data -- normally we 
!!  only evolve on the leaf blocks, but for visualization purposes, 
!!  it is nice to be able to look at other nodes in the tree
!!
!!***

!#define DEBUG_IO

subroutine io_restrictBeforeWrite( restrictNeeded)

#include "constants.h"
  use Driver_interface, ONLY : Driver_getSimTime
  use IO_data, ONLY : io_outputInStack, io_alwaysRestrictCheckpoint
  use Logfile_interface, ONLY : Logfile_stamp

  implicit none
  logical, intent(OUT) :: restrictNeeded

  character (len=500) :: log_message
  real          :: lsimTime   !local sim time variable
  real, save    :: oldsimTime = 0.0  !initial value does not matter
  integer       :: lsimGen    !local sim generation variable
  integer, save :: oldsimGen = -1 !different from any valid generation, to trigger update at startup.


  call Driver_getSimTime(lsimTime, lsimGen)

  !The tree data is already consistent (and so no restriction needed) if:
  !1.  We have already restricted in *this* time step AND
  !2.  The grid has not refined/derefined in *this* time step.
  if (oldsimTime == lsimTime .and. oldsimGen == lsimGen) then
     ![oldSimTime == lSimTime] - same time step.
     ![oldSimGen == lSimGen] - same generation (i.e. no new regrid events).
     restrictNeeded = .false.
  else

     !Set restrictNeeded to .false. if the user has set "io_alwaysRestrictCheckpoint" 
     !to .false. in their flash.par.  This runtime parameter is only 
     !considered when IO_writeCheckpoint and/or IO_writePlotfile are called 
     !directly (i.e. not through IO_Output, IO_OutputInitial or IO_OutputFinal).
     if (.NOT. io_outputInStack .AND. .NOT. io_alwaysRestrictCheckpoint) then
        restrictNeeded = .FALSE.
     else
        !We must restrict.  Save the current simulation time and generation in 
        !static variables.
        restrictNeeded = .true.
        oldsimTime=lsimTime
        oldsimGen =lsimGen
     end if
  end if

#ifdef DEBUG_IO
     write (log_message,'(a,l2,a,l2,a,l2)') "Restrict=", restrictNeeded, &
          ", Called by IO_Output family=", io_outputInStack, &
          ", io_alwaysRestrictCheckpoint=", io_alwaysRestrictCheckpoint
     call Logfile_stamp( log_message, "[io_restrictBeforeWrite]")
#endif

end subroutine io_restrictBeforeWrite
